import os
import torch
import torch.utils.data as data
import numpy as np
from src.pack_model.model import DecoderLSTM, Signal_DataSet, Signal_DataCollection
import random
import multiprocessing
from multiprocessing.pool import Pool
import logging
import traceback
import sys


class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)


class NoDaemonPool(Pool):
    Process = NoDaemonProcess


def predict_one_bin(bin_index, bin_files, device, options):

    try:
        logging.info("Predicting parallel bin index {}".format(bin_index))

        snarl_prediction_res = {}  # a 3-level dict. {snarl_id: {path1: {ref2alt: [], alt2ref: []}, path2: {}}}

        batch_size = 128
        model_architecture = "LSTM"
        model_layer = 1
        model_fc = 64
        model_bidirect = True

        if device == "cpu":
            torch.set_num_threads(options.model_cpu)

        predict_dataset = Signal_DataSet(bin_files, data_type="predict")

        predict_loader = data.DataLoader(predict_dataset, batch_size=batch_size, shuffle=False, collate_fn=Signal_DataCollection, num_workers=2, pin_memory=True)

        model = DecoderLSTM(n_input=predict_dataset.data_dim, n_hidden=model_fc, n_layer=model_layer, architecture=model_architecture, bidirect=model_bidirect).to(device)
        model.load_state_dict(torch.load(options.model_path, map_location=torch.device(device)))
        model.eval()

        for batch_idx, (ids, Xs, Ys, X_lens) in enumerate(predict_loader):
            # # run prediction
            Xs = Xs.to(device)

            pred = model(Xs.to(device), X_lens)

            pred_labels = torch.max(pred, 2)[1].cpu().data.numpy()

            for index in range(len(ids)):
                snarl_prediction_res[ids[index]] = pred_labels[index][: X_lens[index]]

        save_npz_path = os.path.join(options.output_path, "tmp.matrix.{}.res.npz".format(bin_index))
        np.savez(save_npz_path, **snarl_prediction_res)

        return "finish"

    except:

        error_type, error_value, error_trace = sys.exc_info()
        error_log = "Error log: " + str(error_type) + ": " + str(error_value) + '. ' + 'Locate At: ' + str(traceback.extract_tb(error_trace))
        error_info = [bin_index, error_log]

        return error_info


def predict_projection_matrix_parallel(parallel_bin_num, options):

    # # STEP: distribute all the saved npz
    all_npz_files = []
    for file in os.listdir(options.npz_output_path):
        if file.endswith(".npz"):
            all_npz_files.append(os.path.join(options.npz_output_path, file))

    random.shuffle(all_npz_files)

    split_npz_files = [all_npz_files[i::parallel_bin_num] for i in range(parallel_bin_num)]

    # # STEP: run predict
    if options.model_device == "gpu":
        if torch.cuda.is_available():
            device = "cuda"
        else:
            logging.warning("GPU is not available, using CPU instead")
            device = "cpu"
    else:
        device = "cpu"

    errored_bin = []

    if device == "cuda":
        for bin_index in range(parallel_bin_num):
            res = predict_one_bin(bin_index, split_npz_files[bin_index], device, options)

            if res != "finish":
                errored_bin.append(res[0])

    else:
        # for bin_index in range(parallel_bin_num):
        #     predict_one_bin(bin_index, device, options)

        parallel_pool = NoDaemonPool(max(1, int(options.thread_num / options.model_cpu)))
        parallel_pool_res = []

        for bin_index in range(parallel_bin_num):
            # predict_one_bin(bin_index, split_npz_files[bin_index], device, options)
            parallel_pool_res.append(parallel_pool.apply_async(predict_one_bin, (bin_index, split_npz_files[bin_index], device, options)))

        parallel_pool.close()
        parallel_pool.join()

        # STEP: deal with errors from process pool
        for res in parallel_pool_res:
            res = res.get()

            if res != "finish":
                logging.error("Fail predicting for bin {}. {}".format(res[0], res[1]))
                errored_bin.append(res[0])

    if len(errored_bin) != 0:
        logging.error("Exit due to predicting error")
        exit(-1)
#
# def manual_predict():
#     from model import DecoderLSTM, Signal_DataSet, Signal_DataCollection
#
#     # # define the model
#     model_path = "/mnt/d/workspace/svasm/train_sim/new_ssv_csv_k30_s30/LSTM-l2-fc256-bi.pth"
#     data_path = "/mnt/d/workspace/svasm/hprc_hg002/graph_minigraph_chr1/out_test/tmp.matrix.1.npz"
#     # data_path = "/mnt/d/workspace/svasm/train_sim/new_ssv_csv_k30_s30/val_projections.npz"
#
#     batch_size = 128
#     model_architecture = "LSTM"
#     model_layer = 2
#     model_fc = 256
#     model_bidirect = True
#     device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
#
#     predict_dataset = Signal_DataSet(data_path, data_type="predict")
#     predict_loader = data.DataLoader(predict_dataset, batch_size=batch_size, shuffle=False, collate_fn=Signal_DataCollection, num_workers=2, pin_memory=True)
#
#     # device = "cpu"
#     model = DecoderLSTM(n_input=predict_dataset.data_dim, n_hidden=model_fc, n_layer=model_layer, architecture=model_architecture, bidirect=model_bidirect).to(device)
#     model.load_state_dict(torch.load(model_path))
#     model.eval()
#
#     imprecise_cnt = 0
#     for batch_idx, (id, X, y, X_len) in enumerate(predict_loader):
#
#         X, y = X.to(device), y.to(device)
#
#         pred = model(X, X_len)
#         pred_labels = torch.max(pred, 2)[1].cpu().data.numpy()
#
#         target_labels = y.cpu().data.numpy()
#
#         for id_index in range(len(id)):
#             print(id[id_index])
#             print(pred_labels[id_index])
#             print(target_labels[id_index])
#
#             print(len(np.where(pred_labels[id_index] != target_labels[id_index])[0]))
#
#         print()
#         print()


# if __name__ == '__main__':
#     manual_predict()
