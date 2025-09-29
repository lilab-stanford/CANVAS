from __future__ import annotations

import os
from pathlib import Path
from typing import List, Sequence, Tuple, Union, Optional

import numpy as np
import pandas as pd
import torch
from PIL import Image
from torchvision import transforms as T
from timm.data.constants import IMAGENET_INCEPTION_MEAN, IMAGENET_INCEPTION_STD

from .model import CustomModel

__all__ = ["load_model", "predict_images", "predict_folder"]


# ---------- utils ----------

def _resolve_local_weights(path: Optional[str]) -> Optional[str]:
    """
    Resolve a local weights path.
    - If `path` is None, try env `CANVAS_WEIGHTS`.
    - Return a resolved string path if found; otherwise return None.
    - If a path is provided but not found, raise FileNotFoundError.
    """
    w = path or os.getenv("CANVAS_WEIGHTS", "")
    if not w:
        return None
    wp = Path(w).expanduser().resolve()
    if not wp.exists():
        raise FileNotFoundError(f"Weights file not found: {wp}")
    return str(wp)


def _build_transform(img_size: int) -> T.Compose:
    return T.Compose([
        T.Resize((img_size, img_size)),
        T.ToTensor(),
        T.Normalize(mean=IMAGENET_INCEPTION_MEAN, std=IMAGENET_INCEPTION_STD),
    ])


def _load_rgb(path: Union[str, Path]) -> Image.Image:
    with Image.open(path) as im:
        return im.convert("RGB")


def _make_stain_normalize_fn(
    reference_image: Union[str, Path],
    method: str = "Vahadane",
):
    """
    Create a callable that applies stain normalization using tiatoolbox.
    Lazily imported so tiatoolbox is only required when enabled.
    """
    from tiatoolbox.tools import stainnorm

    ref = _load_rgb(reference_image)
    ref_np = np.array(ref)  # HxWx3, uint8
    normalizer = stainnorm.get_normalizer(method)
    normalizer.fit(ref_np)

    def _apply(img_pil: Image.Image) -> Image.Image:
        src = np.array(img_pil)  # uint8 RGB
        dst = normalizer.transform(src)  # uint8 RGB
        return Image.fromarray(dst)

    return _apply


# ---------- public API ----------

def load_model(
    *,
    weights: Optional[str] = None,           # local checkpoint (optional)
    musk_source: str = "hf_hub:xiangjx/musk",
    num_classes: int = 11,
    visual_output_dim: int = 1024,
    device: Optional[str] = None,
    strict: bool = False,
) -> Tuple[torch.nn.Module, str]:
    """
    Build the model and (optionally) load user weights.

    Returns:
        (model, device)

    """
    dev = device or ("cuda" if torch.cuda.is_available() else "cpu")

    # 1) Construct backbone
    model = CustomModel(
        visual_output_dim=visual_output_dim,
        num_classes=num_classes,
        musk_source=musk_source,
    ).to(dev).eval()

    # 2) Optionally load local weights
    wpath = _resolve_local_weights(weights)
    if wpath is None:
        setattr(model, "_canvas_weights_loaded", False)
        return model, dev

    # Safe load
    try:
        state = torch.load(wpath, map_location="cpu", weights_only=True)
    except TypeError:
        state = torch.load(wpath, map_location="cpu")

    loaded = False
    if isinstance(state, dict):
        for key in (None, "state_dict", "model"):
            sd = state if key is None else state.get(key, None)
            if isinstance(sd, dict):
                try:
                    model.load_state_dict(sd, strict=strict)
                    loaded = True
                    break
                except Exception:
                    pass
    if not loaded:
        raise RuntimeError(f"Failed to load weights into model: {wpath}")

    setattr(model, "_canvas_weights_loaded", True)
    return model, dev


def predict_images(
    model: torch.nn.Module,
    device: str,
    image_paths: Sequence[Union[str, Path]],
    *,
    img_size: int = 384,
    batch_size: int = 64,
    color_norm: bool = False,
    reference_image: Optional[Union[str, Path]] = None,
    color_norm_method: str = "Vahadane",
    return_numpy: bool = False,
) -> Union[pd.DataFrame, Tuple[np.ndarray, List[str]]]:
    """
    Pure inference: return per-class probabilities and argmax prediction for each image.

    Optional stain color normalization (default OFF):
      - Set `color_norm=True` AND pass `reference_image=...` to enable.
      - `color_norm_method` defaults to 'Vahadane'.
    """
    tfm = _build_transform(img_size)
    paths = [str(Path(p)) for p in image_paths]
    probs_chunks: List[np.ndarray] = []

    # Optional color normalization
    normalize_fn = None
    if color_norm:
        if not reference_image:
            raise ValueError("color_norm=True requires `reference_image` path.")
        normalize_fn = _make_stain_normalize_fn(reference_image, method=color_norm_method)

    model.eval()
    with torch.no_grad():
        for i in range(0, len(paths), batch_size):
            batch = paths[i:i + batch_size]
            pil_list = [_load_rgb(p) for p in batch]
            if normalize_fn is not None:
                pil_list = [normalize_fn(img) for img in pil_list]
            imgs = [tfm(im) for im in pil_list]
            x = torch.stack(imgs, 0).to(device)
            logits = model(x)  # [N, C]
            prob = torch.softmax(logits, dim=1).cpu().numpy()
            probs_chunks.append(prob)

    probs = np.concatenate(probs_chunks, axis=0)
    if return_numpy:
        return probs, paths

    df = pd.DataFrame(probs, columns=[f"prob_class_{k}" for k in range(probs.shape[1])])
    df.insert(0, "image_path", paths)
    df["pred"] = df[[c for c in df.columns if c.startswith("prob_class_")]].values.argmax(axis=1)
    return df


def predict_folder(
    model: torch.nn.Module,
    device: str,
    folder: Union[str, Path],
    *,
    patterns: Tuple[str, ...] = (".png", ".jpg", ".jpeg", ".tif", ".tiff"),
    recursive: bool = True,
    **kwargs,
):
    """
    Convenience: predict all images under a folder.
    """
    folder = Path(folder)
    if recursive:
        paths = sorted([str(p) for p in folder.rglob("*") if p.suffix.lower() in patterns])
    else:
        paths = sorted([str(p) for p in folder.iterdir() if p.suffix.lower() in patterns])
    return predict_images(model, device, paths, **kwargs)