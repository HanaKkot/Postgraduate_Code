from numpy.core.fromnumeric import shape, size
import torch
import numpy as np
import pickle
import image_processing


from image_processing import read_path, resize_image

import torch.nn.functional as functional
from torch.nn.modules.loss import CrossEntropyLoss


from torchvision import transforms
from torch import nn, optim
from torch.utils.data import Dataset, DataLoader
from PIL import Image

path_name = r'D:\Code1\train'
path_name_hammond = r'D:\Code1\train2'
path_name_ayoan = r'D:\Code1\train3'
path_name_ben = r'D:\Code1\train4'


hana_images, hana_labels = read_path(path_name, 0)
hammond_images, hammond_labels= read_path(path_name_hammond, 1)
ayoan_images, ayoan_labels = read_path(path_name_ayoan, 2)
ben_images, ben_labels = read_path(path_name_ben, 3)

whole_images = hana_images + hammond_images + ayoan_images + ben_images
whole_labels = hana_labels + hammond_labels + ayoan_labels + ben_labels



custom_transforms = transforms.Compose(
    [transforms.ToTensor(),
    transforms.Normalize((0.5,), (0.5,))])

batch_size = 4
classes = ('Hana', 'Ben', 'Hammond', 'Ayoan')

class FacialDatabase(Dataset):
    def __init__(self, features, labels, transform=None):
        self.features = features
        self.labels = labels
        self.transform = transform

    def __len__(self):
        return len(self.features)

    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()

        x = self.features[idx]
        y = self.labels[idx]

        if self.transform:
            x = self.transform(x)
            y = y

        return x, y

trainset = FacialDatabase(whole_images, whole_labels, custom_transforms)
trainloader = DataLoader(trainset, batch_size, shuffle=True)


#create a CNN module
class CNN(nn.Module):
    def __init__(self):
        self.output_size = 4 #Ben, Hana, Hammond, Ayoan
        super(CNN, self).__init__()
        self.conv1 = nn.Conv2d(1, 32, 5)
        self.pool = nn.MaxPool2d(2, 2)
        self.conv2 = nn.Conv2d(32, 64, 5)
        self.conv3 = nn.Conv2d(64, 128, 5)
        self.ful_con1 = nn.Linear(128 * 4 * 4, 120)
        self.ful_con2 = nn.Linear(120, 84)
        self.ful_con3 = nn.Linear(84, self.output_size)
    
    def forward(self, x):
        x = self.pool(functional.relu(self.conv1(x)))
        x = self.pool(functional.relu(self.conv2(x)))
        x = self.pool(functional.relu(self.conv3(x)))
        #print(x.shape)
        x = x.view(-1, 128 * 4 * 4)
        x = functional.relu(self.ful_con1(x))
        x = functional.relu(self.ful_con2(x))
        x = self.ful_con3(x)
        return x
cnn = CNN()




criterion = nn.CrossEntropyLoss()
optimize = optim.Adam(cnn.parameters(), lr = 1e-3)

#train 
'''
for epoch in range(3):

    running_loss = 0.0
    for i, data in enumerate(trainloader, 0):
        inputs, labels = data
        optimize.zero_grad()
        outputs = cnn(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        optimize.step()
        running_loss += loss.item()

        if i % 500 == 499:
            print("Epoch: %d, %5d loss: %.5f" % (epoch + 1, i + 1, running_loss / 500))
            running_loss = 0.0
'''
#torch.save(cnn, "facial_recognition.pt")
#print("Model has been saved sucessfully!")

def cnn_output(image):
    image = resize_image(image)
    image = image.reshape((1 , 1, 64, 64))
    image = image.astype('float32')
    test_dataset = FacialDatabase(image, 'image')
    test_trainloader = DataLoader(test_dataset)
    image /= 255
    model = torch.load("facial_recognition.pt")
    model.eval()    
    
    dataiter = iter(test_trainloader)
    img, lbl = dataiter.next()
    output = model(img)

    sm = nn.Softmax(dim = 1)
    sm_output = sm(output)
    print(sm_output)

    probs, index = torch.max(sm_output, dim =1)
    return probs, index

