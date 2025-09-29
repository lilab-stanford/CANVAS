
import os
from os.path import join
import pandas as pd
import numpy as np
from pathlib import Path
from PIL import Image
import tqdm
from sklearn.metrics import precision_score, recall_score,confusion_matrix

import torch
from torchvision import transforms
from torch.utils.data import Dataset,DataLoader
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import WeightedRandomSampler
from torch.utils.tensorboard import SummaryWriter
from musk import utils
from timm import create_model
from timm.data.constants import IMAGENET_INCEPTION_MEAN, IMAGENET_INCEPTION_STD

from .models import CustomModel

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

def seed_torch(seed=1):
    import random
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if device.type == 'cuda':
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)  # if you are using multi-GPU.
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
seed_torch(1)



class PatchDataset(Dataset):
    def __init__(self, csv, transform=None):
        self.images = csv['image_path'].values
        self.labels = csv['label'].values
        self.transform = transform

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        image_path=self.images[idx]
        label=self.labels[idx]
        image = Image.open(image_path)
        if self.transform:
            image = self.transform(image)
        return image, label,image_path

def print_network(net):
    num_params = 0
    num_params_train = 0
    print(net)

    print("\nTrainable parameters:")
    for name, param in net.named_parameters():
        n = param.numel()
        num_params += n
        if param.requires_grad:
            num_params_train += n
            print(f"{name}, Shape: {param.shape}")

    print('\nTotal number of parameters: %d' % num_params)
    print('Total number of trainable parameters: %d' % num_params_train)


data_dir=""
val_csvs = pd.read_csv(join(data_dir, 'val_labels.csv'))

img_size = 384
transform_val = transforms.Compose([
    transforms.Resize((img_size, img_size)),
    transforms.ToTensor(),
    transforms.Normalize(mean=IMAGENET_INCEPTION_MEAN, std=IMAGENET_INCEPTION_STD)
])

val_dataset = PatchDataset(val_csvs,transform_val)

#data loader
num_workers = 8
val_loader = DataLoader(val_dataset, batch_size=64, shuffle=False, num_workers=num_workers)

num_classes = 11
# Initialize the model
model = CustomModel(visual_output_dim=1024, num_classes=num_classes).to(device)
print_network(model)


check_point_path=f'../checkpoints/checkpoint_epoch_{i}.pth'
checkpoint = torch.load(check_point_path)
model.load_state_dict(checkpoint)
# Validation phase
model.eval()
# Create lists to store image paths and probabilities
all_image_paths = []
all_probs = {f'prob_class_{i}': [] for i in range(11)}
all_labels = []
with torch.no_grad():
    val_loop = tqdm.tqdm(enumerate(val_loader), total=len(val_loader))
    for i, data in val_loop:
        inputs, labels = data[0].to(device), data[1].to(device)
        image_paths = data[2]
        outputs = model(inputs)
        # Get predicted classes
        _, predicted = torch.max(outputs, 1)
        #get probabilities
        probs = torch.nn.functional.softmax(outputs, dim=1)
        all_image_paths.extend(image_paths)
        # Store the probabilities for each class in separate columns
        probs = probs.cpu().numpy()  # Move probs to CPU and convert to numpy array
        for i in range(11):  # Assuming there are 11 classes
            all_probs[f'prob_class_{i}'].extend(probs[:, i])
        all_labels.extend(labels.cpu().numpy())
# Convert the lists to a pandas dataframe
df = pd.DataFrame({
    'image_path': all_image_paths,
    **all_probs,  # Unpack the probabilities dictionary to separate columns
    'label': all_labels
})
df.to_csv("../results.csv", index=False)