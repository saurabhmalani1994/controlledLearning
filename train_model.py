import numpy as np
import matplotlib.pyplot as plt

import torch
from torch.utils.data import DataLoader
from tqdm.auto import tqdm

import utils
import khammash_repro

def main():
    # Generate the data
    phi, L, _ = khammash_repro.datagen(n_traj=35, sp_per_traj=3)
    # phi, L = khammash_repro.datagen(n_traj=1, sp_per_traj=2)
    phi_out_t0, L_out_t0, phi_out_t1 = khammash_repro.time_embed(phi, L, n_embed=7, n_gap=2)

    # Define the dataset
    dataset = khammash_repro.Dataset(phi_out_t0, L_out_t0, phi_out_t1)

    # Train Validation Split
    train_size = int(0.8 * len(dataset))
    val_size = len(dataset) - train_size
    train_dataset, val_dataset = torch.utils.data.random_split(dataset, [train_size, val_size])

    # Define the dataloader
    train_dataloader = DataLoader(train_dataset, batch_size=32, shuffle=True)
    val_dataloader = DataLoader(val_dataset, batch_size=32, shuffle=True)

    # Define the network
    network = utils.MLP(phi_out_t0.shape[-1]+L_out_t0.shape[-1], [100, 100, 100], phi_out_t1.shape[-1])
    network = utils.Network(network, x_norm=[0, 1], p_norm=[0, 800])
    network = network.float()
    if torch.cuda.is_available():
        network = network.cuda()

    # Define the optimizer
    optimizer = torch.optim.Adam(network.parameters(), lr=1e-3)

    # Define the loss function
    loss_func = torch.nn.MSELoss()

    # Define the learning rate scheduler
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, patience=50, factor=0.5, min_lr=0.000001)

    # Total epochs
    epochs = 2000

    # Define the lists to store the training and validation loss
    epoch_arr = []
    epoch_train_loss = []
    epoch_val_loss = []
    epoch_lr = []

    # Define progress bar
    pbar = tqdm(range(0, epochs), leave=True, position=0, desc=progress(0,0))

    # Train the network
    for epoch in pbar:
        epoch_arr.append(epoch)
        cnt, sum_loss = 0, 0
        for x, p, y in train_dataloader:
            if torch.cuda.is_available():
                x = x.float().cuda()
                p = p.float().cuda()
                y = y.float().cuda()

            y_pred = network(x, p)

            loss = loss_func(y_pred, y)
            sum_loss += loss.detach().cpu().numpy()
            cnt += 1

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        epoch_train_loss.append(sum_loss/cnt)
        scheduler.step(sum_loss/cnt)
        epoch_lr.append(scheduler._last_lr)

        cnt, sum_loss = 0, 0
        for x, p, y in val_dataloader:
            if torch.cuda.is_available():
                x = x.float().cuda()
                p = p.float().cuda()
                y = y.float().cuda()

            y_pred = network(x, p)

            loss = loss_func(y_pred, y)
            sum_loss += loss.detach().cpu().numpy()
            cnt += 1

        epoch_val_loss.append(sum_loss/cnt)
        pbar.set_description(progress(epoch_train_loss[-1], epoch_val_loss[-1]))

    # Save the network
    torch.save(network.state_dict(), 'network.net')
    np.save('epoch_arr.npy', epoch_arr)
    np.save('epoch_train_loss.npy', epoch_train_loss)
    np.save('epoch_val_loss.npy', epoch_val_loss)

def progress(train_loss, val_loss):
    """Define progress bar description."""
    return "Train/Loss: {:.2e}  Val/Loss: {:.2e}".format(
        train_loss, val_loss)

if __name__ == "__main__":
    main()