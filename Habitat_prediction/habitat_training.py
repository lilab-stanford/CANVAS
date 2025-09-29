
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

# Custom Model with Additional Layers in the Classification Head
class CustomModel(nn.Module):
    def __init__(self, visual_output_dim, num_classes):
        super(CustomModel, self).__init__()
        model_config = "musk_large_patch16_384"
        model_musk = create_model(model_config, vocab_size=64010)
        checkpoint_path = "/home/usr/.cache/model.safetensors"
        utils.load_model_and_may_interpolate(checkpoint_path, model_musk, 'model|module', '')
        self.visual = model_musk
        self.classification_head = nn.Sequential(
            nn.Linear(visual_output_dim, 256),
            nn.ReLU(),
            nn.Dropout(p=0.5),
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Dropout(p=0.5),
            nn.Linear(128, num_classes)
        )

    def forward(self, x):
        x = self.visual(
            image=x,
            with_head=False,
            out_norm=False
        )[0]
        x = self.classification_head(x)
        return x

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
train_csvs = pd.read_csv(join(data_dir, 'train_labels.csv'))
val_csvs = pd.read_csv(join(data_dir, 'val_labels.csv'))

# Calculate class frequencies
class_counts = train_csvs['label'].value_counts()
class_freq = class_counts / len(train_csvs)  # Frequency of each class
weights = 1.0 / class_freq[train_csvs['label']].values
sampler = WeightedRandomSampler(weights, num_samples=len(weights), replacement=True)

img_size = 384
transform_train = transforms.Compose([
    transforms.RandomHorizontalFlip(),
    transforms.RandomVerticalFlip(),
    transforms.RandomRotation(90),
    transforms.Resize((img_size, img_size)),
    transforms.ColorJitter(brightness=0.1, contrast=0.1, saturation=0.1, hue=0.01),
    transforms.ToTensor(),
    transforms.Normalize(mean=IMAGENET_INCEPTION_MEAN, std=IMAGENET_INCEPTION_STD)
])
transform_val = transforms.Compose([
    transforms.Resize((img_size, img_size)),
    transforms.ToTensor(),
    transforms.Normalize(mean=IMAGENET_INCEPTION_MEAN, std=IMAGENET_INCEPTION_STD)
])

train_dataset = PatchDataset(train_csvs,transform_train)
val_dataset = PatchDataset(val_csvs,transform_val)

#data loader
num_workers = 8
train_loader = DataLoader(train_dataset, batch_size=64, sampler=sampler, num_workers=num_workers)
val_loader = DataLoader(val_dataset, batch_size=64, shuffle=False, num_workers=num_workers)
num_classes = 11

# Initialize the model
model = CustomModel(visual_output_dim=1024, num_classes=num_classes).to(device)
# Freeze the pre-trained layers
for param in model.parameters():
    param.requires_grad = False
# Unfreeze the classification head
for param in model.classification_head.parameters():
    param.requires_grad = True
print_network(model)
# Define the criterion and optimizer
criterion = FocalLoss(alpha=0.25, gamma=2.0, reduction="mean")
optimizer = optim.Adam(model.parameters(), lr=1e-4)
scheduler = torch.optim.lr_scheduler.ExponentialLR(optimizer, gamma=0.95)

# Training Loop
num_epochs = 60
losses = []  # List to store losses
epoch_losses = []  # List to store average losses per epoch
val_losses = []  # List to store validation losses
writer_dir = join("./runs", "habitat_training")
if not os.path.isdir(writer_dir):
    os.mkdir(writer_dir)
writer = SummaryWriter(writer_dir)
checkpoint_dir='../checkpoints/habitat_training'
if not os.path.isdir(checkpoint_dir):
    os.mkdir(checkpoint_dir)

for epoch in range(num_epochs):
    if epoch == 10:
        #train the classification head and backbone
        for layer in model.visual.beit3.encoder.layers[-2:]:
            for param in layer.parameters():
                param.requires_grad = True
        # Unfreeze the final layer norm
        for param in model.visual.beit3.encoder.layer_norm.parameters():
            param.requires_grad = True
    model.train()
    running_loss = 0.0
    all_preds=[]
    all_labels=[]
    train_loop = tqdm.tqdm(enumerate(train_loader), total=len(train_loader))
    for i, data in train_loop:
        inputs, labels = data[0].to(device), data[1].to(device)

        optimizer.zero_grad()
        outputs= model(inputs)
        loss = criterion(outputs, labels)

        # Backward pass and optimization
        loss.backward()
        optimizer.step()

        running_loss += loss.item()

        # Get predicted classes
        _, predicted = torch.max(outputs, 1)

        all_labels.extend(labels.cpu().numpy())
        all_preds.extend(predicted.cpu().numpy())

    # End of epoch: step the LR scheduler
    scheduler.step()
    current_lr = optimizer.param_groups[0]['lr']
    print(f"[Epoch {epoch + 1:02d}] Learning Rate: {current_lr:.2e}")
    # Calculate average loss per epoch and store
    epoch_loss = running_loss / len(train_loader)
    epoch_losses.append((epoch + 1, epoch_loss))
    writer.add_scalar('Loss/train', epoch_loss, epoch + 1)

    # Calculate confusion_matrix
    conf_matrix = confusion_matrix(all_labels, all_preds,labels=range(11))
    # Calculate recall metrics
    tp_per_class = np.diag(conf_matrix)
    all_per_class = np.sum(conf_matrix, axis=1)
    recall_all = tp_per_class / all_per_class

    # Log class-wise training accuracies with TensorBoard
    for j in range(num_classes):
        writer.add_scalar(f'Accuracy/train_class_{j}', recall_all[j], epoch + 1)

    # Log average training accuracy for the epoch
    avg_train_accuracy = np.mean(recall_all)
    writer.add_scalar('Accuracy/train_avg', avg_train_accuracy, epoch + 1)

    # Validation phase
    model.eval()
    val_loss = 0.0
    correct_preds = np.zeros(num_classes)
    total_preds = np.zeros(num_classes)
    all_labels = []
    all_preds = []
    with torch.no_grad():
        val_loop = tqdm.tqdm(enumerate(val_loader), total=len(val_loader))
        for i, data in val_loop:
            inputs, labels = data[0].to(device), data[1].to(device)
            outputs = model(inputs)
            loss = criterion(outputs, labels)
            val_loss += loss.item()
            # Get predicted classes
            _, predicted = torch.max(outputs, 1)

            all_labels.extend(labels.cpu().numpy())
            all_preds.extend(predicted.cpu().numpy())
    # Calculate confusion_matrix
    conf_matrix = confusion_matrix(all_labels, all_preds,labels=range(11))
    # Calculate recall metrics
    tp_per_class = np.diag(conf_matrix)
    all_per_class = np.sum(conf_matrix, axis=1)
    recall_all = tp_per_class / all_per_class

    # Log class-wise validation accuracies with TensorBoard
    for i in range(num_classes):
        writer.add_scalar(f'Accuracy/val_class_{i}', recall_all[i], epoch + 1)

    # Calculate and log average validation loss and accuracy
    val_loss /= len(val_loader)
    val_losses.append((epoch + 1, val_loss))
    writer.add_scalar('Loss/val', val_loss, epoch + 1)

    avg_val_accuracy = np.mean(recall_all)
    writer.add_scalar('Accuracy/val_avg', avg_val_accuracy, epoch + 1)

    torch.save(model.state_dict(), join(checkpoint_dir, f'checkpoint_epoch_{epoch + 1}.pth'))
    print(f"Model weights saved for epoch {epoch + 1}")

print("Finished Training")