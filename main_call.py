import time

from src.pack_graph import op_gfa
from src.pack_dotplot.op_plot import generate_projection_matrix_for_snarls_parallel
from src.pack_model.predict import predict_projection_matrix_parallel
from src.pack_sv.op_sv import interpret_prediction_to_variant_parallel
import logging
import os
import pickle
import shutil


def swave_call(gfa_obj, options):

    if gfa_obj is None:
        gfa_obj = op_gfa.load_pangraph(options)

        # # STEP: load the gfa
        # gfa_pickle_path = os.path.join(options.output_path, "tmp.gfa_minigraph.pkl")
        #
        # if not os.path.exists(gfa_pickle_path):
        #     gfa_obj = op_gfa.load_pangraph(options)
        #     with open(gfa_pickle_path, "wb") as pickle_out:
        #         pickle.dump(gfa_obj, pickle_out)
        #
        # else:
        #     with open(gfa_pickle_path, "rb") as pickle_in:
        #         gfa_obj = pickle.load(pickle_in)

    parallel_bin_num = options.thread_num

    # # STEP: generate
    logging.info("****************** Step1 Generating ******************")
    generate_projection_matrix_for_snarls_parallel(gfa_obj, options)

    # # STEP: predict
    logging.info("****************** Step2 Predicting ******************")
    predict_projection_matrix_parallel(parallel_bin_num, options)

    # # STEP: interpret
    logging.info("****************** Step3 Interpreting ******************")
    interpret_prediction_to_variant_parallel(gfa_obj, parallel_bin_num, options)

    logging.info("****************** Step4 Cleaning ******************")
    logging.info("Clean tmp files")
    shutil.rmtree(options.npz_output_path)
    shutil.rmtree(options.img_output_path)
    os.system("rm {}/tmp.*".format(options.output_path))

    logging.info("****************** Step5 Finishing ******************")

