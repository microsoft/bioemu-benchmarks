
from bioemu.sample import main as sample
import argparse
from bioemu_benchmarks.benchmarks import Benchmark
from bioemu_benchmarks.samples import IndexedSamples, filter_unphysical_samples, find_samples_in_dir
from bioemu_benchmarks.evaluator_utils import evaluator_from_benchmark
import pandas as pd
import os
from glob import glob

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output_dir",
        type=str,
        default="./bioemu_samples",
        help="Directory to save the sampled PDB, XTC files.",    
    )
    parser.add_argument(
        "--test_cases",
        type=str,
        default=".",   
        help="Directory containing the test cases for the benchmarks.",
    )
    parser.add_argument(
        "--prediction_dir",
        type=str,
        default=".",
        help="Path to the CSV file containing predicted sequences.",
    )
    parser.add_argument(
        "--msa_dir",
        type=str,
        default=".",
        help="Path to the directory containing MSA files for each pair. (.a3m) alternative to prediciton dir with sequences",
    )
    parser.add_argument(
        "--model_weights",
        type=str,
        default=".",
        help="Path to the directory with model weights + config used for inference.",
    )
    parser.add_argument(
        "--num_samples",
        type=int,
        default=1,
        help="Number of samples to generate for each sequence.",
    )
    parser.add_argument(
        "--sample_idx",
        type=int,
        default=0,
        help="Index of the sample to use from the MSA files.",
    )
    args = parser.parse_args()
    RESIDUE_LETTERS = [
        "A",
        "R",
        "N",
        "D",
        "C",
        "Q",
        "E",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
        "X",
    ]

    
    """ if args.prediction_dir !=".":
         # load the prediction CSV file
        domain_motion_df = pd.read_csv(os.path.join(args.prediction_dir, "domain_motion_inference.csv"))
        pred_sequences = domain_motion_df['pred_sequence'].to_numpy()
        ids = domain_motion_df['id'].to_numpy()
        for (seq, id) in zip(pred_sequences, ids):
            sub_dir = args.output_dir + "/" + id
            seq = "".join([i if i != 'X' else "A" for i in seq])
            sample(sequence= seq, num_samples=args.num_samples, output_dir=sub_dir, filter_samples= True, batch_size_100 = 20)"""

    # can also pass MSAs
    msa_dir = glob(os.path.join(args.msa_dir, "*.a3m"))

    msa_list = list(msa_dir)
    msa = msa_list[args.sample_idx]

    sub_dir = args.output_dir + "/" + msa.split("/")[-1].split(".")[0]

    sample(sequence= msa, num_samples=args.num_samples, output_dir=sub_dir, filter_samples= True, batch_size_100 = 20, ckpt_path = args.model_weights + "/" + "checkpoint.ckpt", model_config_path= args.model_weights + "/" + "config.yaml")


    # load testcases based on benchmarking set
    testcases = pd.read_csv(os.path.join(args.test_cases, "testcases.csv"))
    testcases = testcases["test_case"].tolist()
    
    benchmark = Benchmark.MULTICONF_DOMAINMOTION

    sequence_samples = find_samples_in_dir(args.output_dir)
    samples = IndexedSamples.from_benchmark(benchmark=benchmark, sequence_samples=sequence_samples, test_case_list=testcases, include_relevant_sequences=False)

    # Filter unphysical-looking samples from getting evaluated
    samples, _sample_stats = filter_unphysical_samples(samples)

    evaluator = evaluator_from_benchmark(benchmark=benchmark)

    results = evaluator(samples)
    results.plot(args.output_dir) 
    results.save_results(args.output_dir)