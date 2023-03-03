import numpy as np
import matplotlib.pyplot as plt

from tqdm.auto import tqdm

import numbers
import math
import torch
import torch.nn as nn

from torch.optim.lr_scheduler import _LRScheduler
torch.set_default_tensor_type(torch.FloatTensor)
torch.set_default_dtype(torch.float)

from torch import Tensor

np.random.seed(1234)

class MLP(nn.Module):
    def __init__(self, input_dim, hidden_cells, output_dim):
        super(MLP, self).__init__()
        
        self.layers = nn.ModuleList()
        self.activation = nn.SiLU()
        
        self.layers.append(nn.Linear(input_dim,hidden_cells[0]))
#         nn.init.kaiming_normal_(self.layers[-1].weight, mode='fan_in', nonlinearity='relu')
#         nn.init.constant_(self.layers[-1].bias, 0)
        for i in range(len(hidden_cells)-1):
            self.layers.append(nn.Linear(hidden_cells[i],hidden_cells[i+1]))
#             nn.init.kaiming_normal_(self.layers[-1].weight, mode='fan_in', nonlinearity='relu')
#             nn.init.constant_(self.layers[-1].bias, 0)
            
        self.output_layer = nn.Linear(hidden_cells[-1], output_dim)
#         nn.init.xavier_normal_(self.output_layer.weight, gain=nn.init.calculate_gain('linear'))
#         nn.init.constant_(self.output_layer.bias, 0)
        
    def forward(self, input):
        for layer in self.layers:
            input = layer(input)
            input = self.activation(input)
        output = self.output_layer(input)
        return output

class Network(nn.Module):
    def __init__(self, network, x_norm, p_norm):
        super(Network, self).__init__()
        
        self.net = network
        self.x_norm_func = lambda x: (x - x_norm[0]) / (x_norm[1] - x_norm[0])
        self.p_norm_func = lambda p: (p - p_norm[0]) / (p_norm[1] - p_norm[0])
        self.x_denorm_func = lambda x: x * (x_norm[1] - x_norm[0]) + x_norm[0]

    def forward(self, x_input, p_input):
        x_input = self.x_norm_func(x_input)
        p_input = self.p_norm_func(p_input)
        input = torch.cat((x_input, p_input), dim=-1)
        output = self.net(input)
        output = self.x_denorm_func(output)
        return output


class Model_Train():
    def __init__(self, dataloader_train, dataloader_val, network, learning_rate=0.01, device=None):
        if torch.cuda.is_available() and device is None:
            self.device = 'cuda'
        elif not torch.cuda.is_available() and device is None:
            self.device = 'cpu'
        else:
            self.device = device

        print('Using:', self.device)
        
        self.net = network.to(self.device)
        self.norm_func = network.norm_func

        self.trainable_parameters = \
            sum(p.numel() for p in self.net.parameters() if p.requires_grad)

        self.dataloader_train = dataloader_train
        self.dataloader_val = dataloader_val
        
        self.lr = learning_rate
        
        self.optimizer = torch.optim.AdamW(self.net.parameters(),
            lr=self.lr, amsgrad=True, weight_decay=0.01)

        self.criterion = torch.nn.MSELoss(reduction='mean').to(self.device)

        self.train_loss = []
        self.val_loss = []
        
        self.lr_epoch = []
        self.lr_track = []

        self.scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            self.optimizer, patience=50, factor=0.5, min_lr=0.000001)
        
    def train(self, epoch):
        """Train model."""
        self.net.train()
        cnt, sum_loss = 0, 0
        iters = len(self.dataloader_train)
 
        for (x_in, p_in, x_out) in self.dataloader_train:
            self.optimizer.zero_grad()

            x_out_hat = self.net(x_in, p_in)
            
            # Calculate loss
            loss = self.criterion(x_out, x_out_hat)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(self.net.parameters(), 0.1)
            self.optimizer.step()
            
            sum_loss += loss.detach().cpu().numpy()
            cnt += 1

            self.scheduler.step(sum_loss / cnt)
            self.lr_epoch.append(epoch)
            self.lr_track.append(self.scheduler._last_lr)
            
        self.train_loss.append(sum_loss/cnt)
        return sum_loss/cnt
    
    def validate(self, epoch):
        """Validate model."""
        self.net.eval()
        cnt, sum_loss = 0, 0
        with torch.no_grad():
            for (x_in, p_in, x_out) in self.dataloader_val:
                x_out_hat = self.net(x_in, p_in)
                
                # Calculate loss
                loss = self.criterion(x_out, x_out_hat)

                sum_loss += loss.detach().cpu().numpy()
                cnt += 1
        self.val_loss.append(sum_loss/cnt)
        return sum_loss/cnt


    def save_network(self, name):
        """Save network weights and training loss history."""
        filename = name +'.net'
        torch.save(self.net.state_dict(), filename)
        np.save(name+'_training_loss.npy', np.array(self.train_loss))
        np.save(name+'_validation_loss.npy', np.array(self.val_loss))
        np.save(name+'_lr_epoch.npy', np.array(self.lr_epoch))
        np.save(name+'_lr_track.npy', np.array(self.lr_track))
        return name

    def load_network(self, name):
        """Load network weights and training loss history."""
        filename = name + '.net'
        self.net.load_state_dict(torch.load(filename))
        self.train_loss = np.load(name+'_training_loss.npy').tolist()
        self.val_loss = np.load(name+'_validation_loss.npy').tolist()

def progress(train_loss, val_loss):
    """Define progress bar description."""
    return "Train/Loss: {:.2e}  Val/Loss: {:.2e}".format(
        train_loss, val_loss)