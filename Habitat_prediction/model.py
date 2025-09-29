
import os
import torch
import torch.nn as nn

from timm import create_model
from musk import utils, modeling
from huggingface_hub import hf_hub_download


def _resolve_musk_ckpt(musk_source: str, filename: str = "model.safetensors") -> str:
    if musk_source.startswith("hf_hub:"):
        repo_id = musk_source.split(":", 1)[1]
        return hf_hub_download(
            repo_id=repo_id,
            filename=filename,
            local_dir_use_symlinks=False,  # safer on some filesystems
        )
    return os.path.expanduser(musk_source)


class CustomModel(nn.Module):
    """
    Inference-only model:
      - Visual encoder: MUSK (loaded via utils.load_model_and_may_interpolate)
      - Classification head: 3-layer MLP as specified
    """
    def __init__(
        self,
        visual_output_dim: int,
        num_classes: int,
        model_config: str = "musk_large_patch16_384",
        musk_source: str = "hf_hub:xiangjx/musk",  # HF Hub by default; can be a local path
    ):
        super().__init__()

        # 1) Create MUSK backbone
        model_musk = create_model(model_config, vocab_size=64010)

        # 2) Load MUSK weights once (HF Hub) and reuse local cache afterwards
        ckpt_path = _resolve_musk_ckpt(musk_source, filename="model.safetensors")
        utils.load_model_and_may_interpolate(
            ckpt_path, model_musk, "model|module", ""
        )

        self.visual = model_musk

        # 3) Classification head
        self.classification_head = nn.Sequential(
            nn.Linear(visual_output_dim, 256),
            nn.ReLU(),
            nn.Dropout(p=0.5),
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Dropout(p=0.5),
            nn.Linear(128, num_classes),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        feat = self.visual(image=x, with_head=False, out_norm=False)[0]
        logits = self.classification_head(feat)
        return logits
