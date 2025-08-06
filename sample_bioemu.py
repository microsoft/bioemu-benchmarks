
from bioemu.sample import main as sample
import argparse
from bioemu_benchmarks.benchmarks import Benchmark
from bioemu_benchmarks.samples import IndexedSamples, filter_unphysical_samples, find_samples_in_dir
from bioemu_benchmarks.evaluator_utils import evaluator_from_benchmark
import pandas as pd
import os

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output_dir",
        type=str,
        default="./bioemu_samples",
        help="Directory to save the sampled PDB, XTC files.",    
    )
    parser.add_argument(
        "--testcases",
        type=str,
        default=".",   
        help="Directory containing the test cases for the benchmarks.",
    )
    parser.add_argument(
        "--prediction_dir",
        type=str,
        default="dynamicmpnn_infer",
        help="Path to the CSV file containing predicted sequences.",
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
    #save the prediciton sequence as a FASTA file
    # model predict a seq and save it as a FASTA file

    domain_motion_df = pd.read_csv(os.path.join(args.prediction_dir, "domain_motion_inference.csv"))
    pred_sequences = domain_motion_df['pred_sequence'].to_numpy()
    ids = domain_motion_df['id'].to_numpy()
    ids = domain_motion_df['id'].to_numpy()

    # sample a sequence -> get pdb, xtc file . give the seqeunce produced by dynamicmpnn for each bnechmark 
    #seq ="GYDGGAAAAA"
    # can also pass a path to a single sequence FASTA file

    ## check this

    for (seq, id) in zip(pred_sequences, ids):
        sub_dir = args.output_dir + "/" + id
        
        seq = "".join([i if i != 'X' else "A" for i in seq])

        sample(sequence= seq, num_samples=100, output_dir=sub_dir, filter_samples= False)

    # load testcases based on benchmarking set
    testcases = pd.read_csv(os.path.join(args.testcases, "testcases.csv"))
    testcases = testcases["test_case"].tolist()
    
    benchmark = Benchmark.MULTICONF_DOMAINMOTION

    sequence_samples = find_samples_in_dir(args.output_dir)
    samples = IndexedSamples.from_benchmark(benchmark=benchmark, sequence_samples=sequence_samples, test_case_list=testcases, include_relevant_sequences=False)

    # Filter unphysical-looking samples from getting evaluated
    #samples, _sample_stats = filter_unphysical_samples(samples)

    evaluator = evaluator_from_benchmark(benchmark=benchmark)

    results = evaluator(samples)
    results.plot(args.output_dir) 
    results.save_results(args.output_dir)