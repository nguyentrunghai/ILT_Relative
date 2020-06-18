
import numpy as np

YANK_LIGANDS = {}

YANK_LIGANDS["1-methylpyrrole.A__AAA"]      = "methylpyrrole"
YANK_LIGANDS["benzene.A__AAA"]              = "benzene"
YANK_LIGANDS["lysozyme.active.A__AAK"]      = "benzofuran"
YANK_LIGANDS["lysozyme.active.A__AAZ"]      = "allyl ethyl sulfide"
YANK_LIGANDS["lysozyme.active.A__ABA"]      = "hexafluorobenzene"
YANK_LIGANDS["lysozyme.active.A__ABC"]      = "indole"
YANK_LIGANDS["lysozyme.active.A__ABG"]      = "m-xylene"
YANK_LIGANDS["lysozyme.active.A__ABJ"]      = "n-hexylbenzene"
YANK_LIGANDS["lysozyme.active.A__ABK"]      = "nitrobenzene"
YANK_LIGANDS["lysozyme.active.A__AAU"]      = "4-ethyltoluene"

YANK_LIGANDS["lysozyme.inactive.A__AA0"]    = "ethanol"
YANK_LIGANDS["lysozyme.inactive.A__AA1"]    = "methanol"
YANK_LIGANDS["lysozyme.inactive.A__AA3"]    = "dimethyl-sulfoxide"
YANK_LIGANDS["lysozyme.inactive.A__AAS"]    = "DL-camphor"
YANK_LIGANDS["lysozyme.inactive.A__AAY"]    = "1-propanol"
YANK_LIGANDS["lysozyme.inactive.A__AAZ"]    = "1,1-diethylurea"
YANK_LIGANDS["lysozyme.inactive.A__ABA"]    = "1,4 diiodobenzene"
YANK_LIGANDS["lysozyme.inactive.A__ABB"]    = "1,2-diiodobenzene"
YANK_LIGANDS["lysozyme.inactive.A__ABH"]    = "2-bromoethanol"
YANK_LIGANDS["lysozyme.inactive.A__ABI"]    = "2-iodoethanol"
YANK_LIGANDS["lysozyme.inactive.A__ABJ"]    = "benzyl alcohol"
YANK_LIGANDS["lysozyme.inactive.A__ABK"]    = "benzaldehyde oxime"
YANK_LIGANDS["phenol.A__AAA"]               = "phenol"

YANK_LIGANDS["p-xylene.A__AAA"]             = "p-xylene"


def load_scores(file_name, id_col, score_col, std_col, exclude_ligands):
    scores = {}
    standard_devs = {}

    with open(file_name, "r") as handle:
        for line in handle:
            if not line.strip().startswith("#"):
                entries = line.split()
                id = entries[id_col]
                val = entries[score_col]
                std = entries[std_col]

                if id not in exclude_ligands:
                    if val.lower() not in ["inf", "nan"]:
                        scores[id] = np.float(val)
                        standard_devs[id] = np.float(std)
    return scores, standard_devs


def load_exper_bfes(file_name, id_col, score_col, exclude_ligands):
    scores = {}
    with open(file_name, "r") as handle:
        for line in handle:
            if not line.strip().startswith("#"):
                entries = line.split()
                id = entries[id_col]
                val = entries[score_col]

                if id not in exclude_ligands:
                    if val.lower() not in ["inf", "nan"]:
                        scores[id] = np.float(val)
    return scores


def matching_scores(score1, score2, err1, err2, allowed_ligands=[]):
    ligands = set(score1.keys()).intersection(set(score2.keys()))

    if len(allowed_ligands) > 0:
        ligands = ligands.intersection(allowed_ligands)
    ligands = list(ligands)

    x, y, xerr, yerr = [], [], [], []

    for ligand in ligands:
        x.append(score1[ligand])
        y.append(score2[ligand])

        xerr.append(err1[ligand])
        yerr.append(err2[ligand])

    return np.array(x), np.array(y), np.array(xerr), np.array(yerr), ligands


def write_pairs(score1, score2, err1, err2, out, allowed_ligands):
    ligands = set(score1.keys()).intersection( set(score2.keys()))
    if len(allowed_ligands) > 0:
        ligands = ligands.intersection( set(allowed_ligands))
    ligands = list(ligands)

    diff = {l:np.abs(score1[l] - score2[l]) for l in ligands}
    ligands.sort(key=lambda l: diff[l], reverse=True)

    with open(out, "w") as handle:
        handle.write("# ligand     x     y    xerr    yerr   diff\n")
        for l in ligands:
            handle.write("%30s %20.10f %20.10f %20.10f %20.10f %20.10f\n" % (l, score1[l], score2[l],
                                                                             err1[l], err2[l], diff[l]))
    return None

