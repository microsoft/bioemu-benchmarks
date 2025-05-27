import contextlib
import hashlib
import os
from itertools import combinations

import joblib
import mdtraj
import mdtraj.utils
import numpy as np
import requests
from mdtraj import Trajectory
from tqdm import tqdm

StrPath = str | os.PathLike

BIOEMU_HEADER = """
                                  -%=.                           -*                    
                              :@:  ==:::-+%                   :%+                      
                         .=:.  **+   *::::::%              ++*#                        
                          -@+:  ++%:  =:::+%::+          %#+-*      :-                 
                           -*#-  ++++   +:::*=:*       **+=.-   =%=#                   
                   -:==:    :*-% .++=+   #:::==:*     @==:.. .#+=@-::        .*+--+**  
                     %+==.   :#:+--#==-::=#:::+=:+   @==: -+=:=%+.       =+-=+%*       
                      :%==:   :#::+=*=-.+==@:::#::+ @==: @-::#*-     --%@#=**          
                        -#==  .:*::*+-=-:+==*:::+:-@===%+-:-%=.   -*@%+=-+=            
                          #*=: .:+:=*=-=-=*%=-::#:@===%-::-#-  =*#@====+*              
                           +#=- .==:+*=-=:=@==::**#==+=::-=:.%#-%=====#:      :+*+=    
                            +=+=:-*=-+==--:=*=-:*#=+++=:-=+@=:*=-===*=:   =%*:::++     
                 -=-:        %=*=-%-:+-===:-%==:*%=++=:-=@-::%-====%: :*+===+#-        
                   %=+%%.    :#=*=#=---===:-+%==*#=++-:%@:::#===-=%*=:====@.           
                     +*===*+ .%:-+*+==::==-:=%+=*#+++:=*::::+*--+%:-==-=#              
                       :%====%.::-+++-:::==:=+#=*#+++=+=:::==::#=-===:@                
                          ++==*-:::=++:::==:-+@==#=++++-::==:::====:#: .+##*+++#       
                            +==++-:-++=::==-:+@==@=++++::==::::==--@#+-====::=#        
                   ===--::-=-.-===-.:=-:::=-:+%:-#*=+++===::::==:+=:-===::-#.          
                *+:-:--=*%+====-#===:::::::::##:::%=*+++==::===-:-===-:.==             
                  -**#%%##+=*+:-=====:-----:-@=:..=%=++=+=:====:===-:.+:               
                           .+%%#==++=---:-::++:==--%%*#+==-=--+==-:: %:::              
                            =-=*@%@%*%@%+--=*:#-----+@@*=+@#+=+#@%%%%+===:.            
                        .==*+#%+:-+*+-.:--%#:%::::----*=**==-::+%%%@%%%#=:+-           
                      :#*+:=%::%      -*....%::..:::--:@.....*%=:=%*-##*+*=+           
                      :*+=+-#%   :==:.  %..%::...:-----@=..#=.:.    #+--+-*+**         
                       *+**-+   --@@@@=::#%::....---.:--:.@-+@@@%=   %.:**=+#=*        
                        =+#*%   -#@@@@@-:%%::....--:.:--@#++@@@@@%:  %::=-::*- :       
                         =+:=   -*@@@%@:%*--:::.:----==-%==#@@@@@%   @.-:--+#+.        
                           .*#  .-#@@*:+@*::%-:..-==#*--@%@=*%%*-   @:*#-              
                             :+   .:. %%:#@:=:::.---::-+-%@@-      @%.                 
                         =*++*+@%+==#**#-@@%:.:::--::+@@@-==*=#*+#.                    
                        :::::=#-#**+-=@::@-......:-:.::@@:-::%-.--%=++%=.              
                         .-++--::%%--:::::.....::::-:::..:::::-#*%:::-++=#:            
                         *#=--....+@=@+.:.....:--------:..:-@@@#....+=-*:              
                        .--:-=-:-::..#%%-:....:::----=---+@@*:..::-*-::+-              
                        ::::.  .-#@=-..%+#:.:::.:::--::-@%=...%#=::%::+-::             
                               +=-::*=:-*=@-:--:-::::-@=@::+:..:==+:-=:                
                                  *===+++:@##:------@+@:---::+*.                       
                                     :+#=..=#*#::=@%=.:--=%==++++                      
                                     .+%=+-=-.#@@+===--@-:::*                          
                                      :=@+-%====+=++**%##%-@+.                         
                                       @@=-+=%+%+*#*--+%==@#                           
                                       -@@+-:======-=====@@.                           
                                        *@+=-::::-======@@%                            
                                         @+--:::::=====%@@                             
______ _       _____                          ______                 _                          _    
| ___ (_)     |  ___|                         | ___ \               | |                        | |   
| |_/ /_  ___ | |__ _ __ ___  _   _   ______  | |_/ / ___ _ __   ___| |__  _ __ ___   __ _ _ __| | __
| ___ \ |/ _ \|  __| '_ ` _ \| | | | |______| | ___ \/ _ \ '_ \ / __| '_ \| '_ ` _ \ / _` | '__| |/ /
| |_/ / | (_) | |__| | | | | | |_| |          | |_/ /  __/ | | | (__| | | | | | | | | (_| | |  |   < 
\____/|_|\___/\____/_| |_| |_|\__,_|          \____/ \___|_| |_|\___|_| |_|_| |_| |_|\__,_|_|  |_|\_\\
"""


def shahexencode(s: str) -> str:
    """Simple sha256 string encoding"""
    return hashlib.sha256(s.encode()).hexdigest()


def download(src: str, dest: StrPath) -> None:
    """Simple request.get call with a progress bar
    Args:
        src: URL to be retrieved
        dest: Local destination path
    """
    req = requests.get(src, stream=True)
    assert req.status_code == 200
    tsize = int(req.headers.get("content-length", 0))
    progress = tqdm(total=tsize, unit="iB", unit_scale=True, position=0, leave=False)

    with open(dest, "wb") as handle:
        progress.set_description(os.path.basename(dest))
        for chunk in req.iter_content(chunk_size=1024):
            handle.write(chunk)
            progress.update(len(chunk))


@contextlib.contextmanager
def tqdm_joblib(tqdm_object: tqdm):
    """Context manager to patch joblib to report into tqdm progress bar given as argument
    Taken from
    https://stackoverflow.com/questions/24983493/tracking-progress-of-joblib-parallel-execution/58936697#58936697
    """

    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()


def setup_external_dependencies():
    # Some dependencies (e.g., TM-align) need to be setup in advance they're called
    # under the `evaluate_multiconf` function, as systems are evaluated in parallel
    # and the setup can cause race conditions.
    from bioemu_benchmarks.eval.multiconf.align import setup_tm_align

    setup_tm_align()


def filter_unphysical_traj_masks(
    traj: Trajectory,
    max_ca_seq_distance: float = 4.5,
    max_cn_seq_distance: float = 2.0,
    clash_distance: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    See `filter_unphysical_traj` for more details.
    """
    # CA-CA residue distance between sequential neighbouring pairs.
    seq_contiguous_resid_pairs = np.array(
        [(r.index, r.index + 1) for r in list(traj.topology.residues)[:-1]]
    )

    ca_seq_distances, _ = mdtraj.compute_contacts(
        traj, scheme="ca", contacts=seq_contiguous_resid_pairs, periodic=False
    )
    ca_seq_distances = mdtraj.utils.in_units_of(ca_seq_distances, "nanometers", "angstrom")

    frames_match_ca_seq_distance = np.all(ca_seq_distances < max_ca_seq_distance, axis=1)

    # C-N distance between sequential neighbouring pairs.
    cn_atom_pair_indices: list[tuple[int, int]] = []

    for resid_i, resid_j in seq_contiguous_resid_pairs:
        residue_i, residue_j = (
            traj.topology.residue(resid_i),
            traj.topology.residue(resid_j),
        )
        c_i, n_j = (
            list(residue_i.atoms_by_name("C")),
            list(residue_j.atoms_by_name("N")),
        )
        assert len(c_i) == len(n_j) == 1
        cn_atom_pair_indices.append((c_i[0].index, n_j[0].index))

    assert cn_atom_pair_indices

    cn_seq_distances = mdtraj.compute_distances(traj, cn_atom_pair_indices, periodic=False)
    cn_seq_distances = mdtraj.utils.in_units_of(cn_seq_distances, "nanometers", "angstrom")

    frames_match_cn_seq_distance = np.all(cn_seq_distances < max_cn_seq_distance, axis=1)

    # Clashes between any two atoms from different residues
    # Maximum cutoff we need to look at is maximal distance of a frame atom from Ca (O with Ca-O
    # ~ 2.4 Angstrom, 2.5 with small buffer for safety) times two plus the clash limit we set.
    maximal_ca_with_clash = 2 * 2.5 + clash_distance
    rest_distances, _ = _compute_filtered_contacts(
        traj, ca_clash_cutoff=maximal_ca_with_clash, sequence_separation=2
    )
    frames_non_clash = np.all(
        mdtraj.utils.in_units_of(rest_distances, "nanometers", "angstrom") > clash_distance,
        axis=1,
    )
    return frames_match_ca_seq_distance, frames_match_cn_seq_distance, frames_non_clash


def _compute_filtered_contacts(
    trajectory: mdtraj.Trajectory,
    ca_clash_cutoff: float,
    sequence_separation: int = 2,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute close contacts between atoms in residues. Compared to `mdtraj.compute_contacts`, this
    routine first filters residues based on their Ca distances. Only residues with their Ca within
    a cutoff will be considered for the final contact computation, reducing the number of contacts
    considered significantly. This should speed up contact computation by a factor of ten.

    Args:
        trajectory: Trajectory to compute contacts for.
        ca_clash_cutoff: Cutoff used for Ca based prefiltering. Units are Angstrom.
        sequence_separation: Only consider contacts between atoms at least this many residues apart
          in the sequence.

    Returns:
        Contact distances [n_frames x contacts] and indices of pairs [n_contacts x 2]. These will be
        different from the `mdtraj.compute_contacts` due to the filtering, but will lead to the same
        masks (if `sequence_separation` is set to 2).
    """
    # Get indices of Ca atoms and build an array for mapping Ca atom indices to their respective
    # residue indices.
    idx_ca: list[int] = []
    map_ca_to_resid: np.ndarray = np.full((trajectory.top.n_atoms,), -1)
    for atom in trajectory.top.atoms:
        if atom.name == "CA":
            idx_ca.append(atom.index)
            map_ca_to_resid[atom.index] = atom.residue.index

    # Generate all possible Ca index pairs (`combinations` excludes self pairs) and compute Ca
    # distances for preliminary filtering.
    idx_ca_pairs: np.ndarray = np.array(list(combinations(idx_ca, 2)))
    distances_ca = mdtraj.compute_distances(trajectory, atom_pairs=idx_ca_pairs, periodic=False)

    # Convert to Angstrom and filter. The Ca distances still have the trajectory dimension
    # (n_frames x n_distances), for the filter we consider any contact that is closer than the
    # cutoff at one point in the trajectory.
    within_ca_cutoff = np.any(
        mdtraj.utils.in_units_of(distances_ca, "nanometers", "angstrom") < ca_clash_cutoff,
        axis=0,
    )

    # Select contacts within cutoff and convert atom indices to residue indices for contact
    # computation.
    idx_contact_residues = map_ca_to_resid[idx_ca_pairs[within_ca_cutoff]]

    # Native `compute_contacts` only considers residues at least two amino acids apart in the
    # sequence for contact computation. We apply a similar filter here (`sequence_separation=2`)
    # will reproduce the `compute_contacts` behavior.
    mask_sequence_separation = (
        np.abs(idx_contact_residues[:, 0] - idx_contact_residues[:, 1]) > sequence_separation
    )
    idx_contact_residues = idx_contact_residues[mask_sequence_separation]

    # Compute contacts on the filtered residue pairs.
    contact_distances, contact_idx = mdtraj.compute_contacts(
        trajectory, contacts=idx_contact_residues, periodic=False
    )
    return contact_distances, contact_idx


def get_physical_traj_indices(
    traj: Trajectory,
    max_ca_seq_distance: float = 4.5,
    max_cn_seq_distance: float = 2.0,
    clash_distance: float = 1.0,
    strict: bool = False,
) -> np.ndarray:
    """
    See `filter_unphysical_traj`. This returns trajectory frame indices satisfying certain physical criteria.
    """
    (
        frames_match_ca_seq_distance,
        frames_match_cn_seq_distance,
        frames_non_clash,
    ) = filter_unphysical_traj_masks(traj, max_ca_seq_distance, max_cn_seq_distance, clash_distance)
    matches_all = frames_match_ca_seq_distance & frames_match_cn_seq_distance & frames_non_clash
    if strict:
        assert matches_all.sum() > 0, "Ended up with empty trajectory"
    return np.where(matches_all)[0]


def filter_unphysical_traj(
    traj: Trajectory,
    max_ca_seq_distance: float = 4.5,
    max_cn_seq_distance: float = 2.0,
    clash_distance: float = 1.0,
    strict: bool = False,
) -> Trajectory:
    """
    Filters out 'unphysical' frames from a samples trajectory

    Args:
        traj: A trajectory object with multiple frames
        max_ca_seq_distance: Maximum carbon alpha distance between any two contiguous residues in the sequence (in Angstrom)
        max_cn_seq_distance: Maximum carbon-nitrogen distance between any two contiguous residues in the sequence (in Angstrom)
        clash_distance: Minimum distance between any two atoms belonging to different residues (in Angstrom)
        strict: Raises an error if all frames in `traj` are filtered out
    """
    matches_all = get_physical_traj_indices(
        traj=traj,
        max_ca_seq_distance=max_ca_seq_distance,
        max_cn_seq_distance=max_cn_seq_distance,
        clash_distance=clash_distance,
        strict=strict,
    )
    return traj.slice(matches_all, copy=True)
