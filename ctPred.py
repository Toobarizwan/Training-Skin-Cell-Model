import torch
from torch.utils import data
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from sklearn.preprocessing import StandardScaler
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

# specify the device that you'll use cpu or gpu
device = torch.device('cpu')

#1: Pre-process the epigenomics file and expression file 
#Use script

#2: Load, split, and normalize the data
# load the data (format in .csv as an example)
gen_data_p = '/home/tooba/Skin_Cell_Model/processed_IbG/processed_IbG_pseudobulk_Prolif__KC.csv'
gen_data = pd.read_csv(gen_data_p, index_col=0)
gen_data

#Fill NAs with 0 
gen_data['percentile'] = gen_data['percentile'].fillna(0)

# define the training, validation, and test set by chromosomes 
train_set = ["chr1", "chr10", "chr13", "chr15", "chr16", "chr17", "chr18", "chr19", "chr2", "chr21", "chr22", "chr3", "chr4", "chr6", "chr8", "chr9", "chrX", "chrY"]
val_set = ["chr11", "chr14", "chr7"]
test_set = ["chr12", "chr20", "chr5"]


def data_prepare(gen_data, train_set, val_set, test_set, is_normalization=True):
    # split the data
    train_data = gen_data[gen_data['chromo'].isin(train_set)]
    val_data = gen_data[gen_data['chromo'].isin(val_set)]
    test_data = gen_data[gen_data['chromo'].isin(test_set)]


    # choose naive B cell as an example to train a model for gene expression prediction
    train_epi = torch.tensor(train_data.iloc[:, 2:5315].values, dtype=torch.float32).to(device)
    train_exp = torch.tensor(train_data.iloc[:, -1].values, dtype=torch.float32).to(device)

    val_epi = torch.tensor(val_data.iloc[:, 2:5315].values, dtype=torch.float32).to(device)
    val_exp = torch.tensor(val_data.iloc[:, -1].values, dtype=torch.float32).to(device)

    test_epi = torch.tensor(test_data.iloc[:, 2:5315].values, dtype=torch.float32).to(device)
    test_exp = torch.tensor(test_data.iloc[:, -1].values, dtype=torch.float32).to(device)

    if is_normalization:
        # normalize the input data (this step is not necessary but recommanded)
        scaler = StandardScaler().fit(train_epi)
        train_epi = torch.tensor(scaler.transform(train_epi), dtype=torch.float32).to(device)
        val_epi = torch.tensor(scaler.transform(val_epi), dtype=torch.float32).to(device)
        test_epi = torch.tensor(scaler.transform(test_epi), dtype=torch.float32).to(device)

    return train_epi, train_exp, val_epi, val_exp, test_epi, test_exp


# put the data into a iterator which reads data in batches
def dataloader(features, labels, batch_size, is_train = True):
    dataset = data.TensorDataset(features, labels)
    return data.DataLoader(dataset, batch_size, shuffle = is_train)


train_epi, train_exp, val_epi, val_exp, test_epi, test_exp = data_prepare(gen_data, train_set, val_set, test_set)
train_data_iter = dataloader(train_epi, train_exp, batch_size=1000)

len(train_epi)

#3: Model Training and Saving

# define the scPred model and some helper functions

class scPred(nn.Module):
    def __init__(self, **kwargs):
        super().__init__()
        
        scPred_defaults = {
            'num_layers' : 4,
            'input_dim' : 5313,
            'hidden_dim' : 64,
            'output_dim' : 1,
            'reg_lambda' : 5e-4,
            'dropout_rate' : 0.05,
            'learning_rate' : 9e-5,
            'random_seed' : 1024
        }

        scPred_defaults.update(kwargs)

        for key, value in scPred_defaults.items():
            setattr(self, key, value)


        torch.manual_seed(self.random_seed)


        # model main
        layers = [nn.Linear(self.input_dim, self.hidden_dim), nn.ReLU(), nn.Dropout(self.dropout_rate)]
        hidden_layer = [nn.Linear(self.hidden_dim, self.hidden_dim), nn.ReLU(), nn.Dropout(self.dropout_rate)]
        
        for _ in range(self.num_layers - 1):
            layers.extend(hidden_layer)
        
        layers.append(nn.Linear(self.hidden_dim, self.output_dim))
        
        self.net = nn.Sequential(*layers)

    
    def custom_loss(self, y_true, y_pred):
        return F.mse_loss(y_true.reshape(-1, 1), y_pred.reshape(-1, 1))

    def forward(self, x):
        return self.net(x)
    
    def compile(self):
        self.optimizer = optim.Adam(self.parameters(), lr = self.learning_rate, weight_decay = self.reg_lambda)

def plot_loss_curve(train_losses, val_losses, epochs):

    plt.figure(figsize=(8, 6))
    plt.plot(range(1, epochs + 1), train_losses, label='Train_loss')
    plt.plot(range(1, epochs + 1), val_losses, label='Val_loss')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.title('Loss curve')
    plt.legend()
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.show()


def save_model(model, path):
    torch.save(model.state_dict(), path)


def scPred_training(model, train_data_iter, val_epi, val_exp, epochs=80, plot_loss=True, model_path='scPred.pt'):
    model.compile()
    optimizer = model.optimizer

    train_losses = []
    val_losses = []
    best_val_loss = float('inf')
    
    for epoch in range(epochs):
        # Validation
        model.eval()
        with torch.no_grad():
            val_loss = model.custom_loss(model(val_epi), val_exp).item()
        val_losses.append(val_loss)

        # Training
        model.train()
        for batch_x, batch_y in train_data_iter:
            batch_x, batch_y = batch_x.to(device), batch_y.to(device)
            optimizer.zero_grad()
            outputs = model(batch_x)
            loss = model.custom_loss(outputs, batch_y)
            loss.backward()
            optimizer.step()
        
        train_loss_epoch = model.custom_loss(model(train_epi), train_exp).item()
        train_losses.append(train_loss_epoch)

        # Conditional save
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            save_model(model, model_path)
            print(f"Model checkpoint saved at epoch {epoch} with val_loss={val_loss:.4f}")

    # Force final save
    save_model(model, model_path)
    print(f"Final model saved to {model_path}")

    if plot_loss:
        plot_loss_curve(train_losses, val_losses, epochs)

scPred_test = scPred().to(device)
saved_path = '/home/tooba/Skin_Cell_Model/models/scPred_Prolif_KC.pt'
scPred_training(scPred_test, train_data_iter, val_epi, val_exp, epochs=100, model_path=saved_path)


#Evaluate the scPred performance on the test set

def plot_prediction(model, test_features, test_labels):
    png_path = '/home/tooba/Skin_Cell_Model/models/scPred_Prolif_KC.png'
    model.eval()

    plt.figure(figsize=(8, 6))
    plt.scatter(test_labels, model(test_features).data, s=5)
    plt.xlabel('Observed expressions')
    plt.ylabel('Predicted expressions')
    plt.title('Model performance on the test set')

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.show()

# Save the plot
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {png_path}")
    
    # Optionally still show the plot
    plt.show()

plot_prediction(scPred_test, test_epi, test_exp)

#Re-load the training model
def load_model(filepath, model_class):

    model = model_class()
    model.load_state_dict(torch.load(filepath, weights_only = True, map_location = device))

    return model

scPred_reload = load_model('/home/tooba/Skin_Cell_Model/models/scPred_Prolif_KC.pt', scPred)
