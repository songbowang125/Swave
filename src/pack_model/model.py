import numpy as np
import torch
import torch.nn as nn
import torch.utils.data as data

label2index = {"None": 0, "REF": 1, "DEL": 2, "INV": 3, "DUP": 4, "invDUP": 5}
index2label = {0: "None", 1: "REF", 2: "DEL", 3: "INV", 4: "DUP", 5: "invDUP"}

# label2index = {"None": 0, "REF": 1, "DEL": 2}


class Signal_DataSet(data.Dataset):

    def __init__(self, npz_files, data_type):
        self.data_type = data_type

        self.npz_data = []
        for npz_file in npz_files:
            # self.npz_data.extend(np.load(npz_file).items())

            try:
                self.npz_data.extend(np.load(npz_file).items())
            except:
                continue

        self.npz_data = sorted(self.npz_data, key=lambda item: len(item[1]))

        self.npz_data = {k: v for k, v in self.npz_data }

        self.data_matrixs = list(self.npz_data.values())
        self.data_ids = list(self.npz_data.keys())

        # self.data_dim = len(self.data_matrixs[0][0]) - 2
        self.data_dim = 4

    def __getitem__(self, idx):

        data = self.data_matrixs[idx]

        id = self.data_ids[idx]

        seq = torch.as_tensor([[float(seg[1]) - float(seg[0]), float(seg[2]), float(seg[3]), float(seg[4])] for seg in data])

        label = torch.as_tensor([label2index[seg[5]] for seg in data])

        return id, seq, label

    def __len__(self):
        return len(self.data_matrixs)


def Signal_DataCollection(batch):

    # batch.sort(key=lambda x: len(x), reverse=True)

    # # extract sequences and labels
    ids = [seq_label_pair[0] for seq_label_pair in batch]
    sequences = [seq_label_pair[1] for seq_label_pair in batch]
    labels = [seq_label_pair[2] for seq_label_pair in batch]
    lengths = [len(seq) for seq in sequences]

    # # padding
    padded_sequences = nn.utils.rnn.pad_sequence(sequences, batch_first=True)
    padded_labels = nn.utils.rnn.pad_sequence(labels, batch_first=True)

    return ids, padded_sequences, padded_labels, lengths


class DecoderLSTM(nn.Module):
    def __init__(self, n_input=3, n_layer=2, n_hidden=256, n_class=len(list(label2index.keys())), architecture="LSTM", bidirect=True):
        super(DecoderLSTM, self).__init__()

        self.n_input = n_input
        self.n_layer = n_layer  # RNN hidden layers
        self.n_hidden = n_hidden  # RNN hidden nodes
        self.n_class = n_class

        self.architecture = architecture
        self.bidirectional = bidirect

        if self.architecture == "LSTM":
            self.rnn = nn.LSTM(
                input_size=self.n_input, hidden_size=self.n_hidden, num_layers=self.n_layer, batch_first=True, bidirectional=self.bidirectional,
            )

        elif self.architecture == "GRU":
            self.rnn = nn.GRU(
                input_size=self.n_input, hidden_size=self.n_hidden, num_layers=self.n_layer, batch_first=True, bidirectional=self.bidirectional,
            )
        else:
            print("No this architecture", self.architecture)

        #
        # self.fc1 = nn.Linear(self.n_hidden, self.n_fc)
        # self.fc2 = nn.Linear(self.n_fc, self.n_class)

        if self.bidirectional:
            self.fc = nn.Linear(self.n_hidden * 2, self.n_class)
        else:
            self.fc = nn.Linear(self.n_hidden, self.n_class)

    def forward(self, X, X_len):

        self.rnn.flatten_parameters()

        # lstm_out, (h_n, h_c) = self.LSTM(X, None)

        rnn_out, (h_n, h_c) = self.rnn(nn.utils.rnn.pack_padded_sequence(X, lengths=X_len, batch_first=True, enforce_sorted=False), None)
        rnn_out, _ = nn.utils.rnn.pad_packed_sequence(rnn_out, batch_first=True)

        # FC layers
        fc_out = self.fc(rnn_out)
        # fc_out = self.fc(lstm_out[:, -1, :])    # choose RNN_out at the last time step

        # # two fc layers
        # fc_out = self.fc1(lstm_out)
        # fc_out = F.relu(fc_out)
        # fc_out = F.dropout(fc_out, p=self.p_drop, training=self.training)
        # fc_out = self.fc2(fc_out)

        return fc_out
