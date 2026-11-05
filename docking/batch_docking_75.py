#!/usr/bin/env python3
"""
batch_docking_75.py — Full 75-ligand docking against 4GIZ + 4XR8
Uses validated fix_docking_pipeline.py infrastructure.
"""

import subprocess, json, csv, time
from pathlib import Path
from datetime import datetime

WORKSPACE = Path(r"F:\HPV_structure")
LIGAND_CSV = WORKSPACE / "ligands" / "01_candidate_ligand_library.csv"
OUTPUT_DIR = WORKSPACE / "docking_results_full"
OUTPUT_DIR.mkdir(exist_ok=True)

WSL_DISTRO = "Ubuntu-22.04"

RECEPTORS = {
    "4GIZ": {
        "pdbqt": "/mnt/f/HPV_structure/protein_structures/4GIZ_cleaned.pdbqt",
        "center": [-13.26, 65.03, 30.59],
        "size": [88.31, 60.82, 61.46]
    },
    "4XR8": {
        "pdbqt": "/mnt/f/HPV_structure/protein_structures/4XR8_cleaned.pdbqt",
        "center": [-54.1, -4.1, -14.9],
        "size": [55, 47, 51]
    }
}

def win_to_wsl(p):
    p = Path(p)
    return f"/mnt/{p.drive[0].lower()}{str(p).replace(chr(92), '/')[2:]}"

def run_wsl(cmd, timeout=600):
    full = f'wsl -d {WSL_DISTRO} -- bash -c "{cmd}"'
    r = subprocess.run(full, shell=True, capture_output=True, text=True, timeout=timeout)
    return r.returncode, r.stdout, r.stderr

def parse_vina(stdout):
    energies = []
    in_table = False
    for line in stdout.split('\n'):
        if in_table:
            parts = line.split()
            if len(parts) >= 2 and parts[0].isdigit():
                try: energies.append(float(parts[1]))
                except: continue
        if 'mode' in line.lower() and 'affinity' in line.lower():
            in_table = True
    return energies

def dock_one(receptor_name, cfg, lig_id, smiles, name):
    d = OUTPUT_DIR / receptor_name
    d.mkdir(exist_ok=True)
    smi = d / f"{lig_id}.smi"
    pdb = d / f"{lig_id}.pdb"
    pdbqt = d / f"{lig_id}.pdbqt"
    out = d / f"{lig_id}_out.pdbqt"

    smi.write_text(f"{smiles}\t{name}\n")

    # 3D generation
    rc, so, se = run_wsl(
        f"obabel -ismi {win_to_wsl(smi)} -opdb -O {win_to_wsl(pdb)} --gen3d --minimize --ff MMFF94", 120)
    if rc != 0: return {"ligand_id": lig_id, "receptor": receptor_name, "success": False, "error": f"3D: {se[:80]}"}

    # PDBQT conversion
    rc, so, se = run_wsl(
        f"obabel -ipdb {win_to_wsl(pdb)} -opdbqt -O {win_to_wsl(pdbqt)} -p 7.4", 120)
    if rc != 0: return {"ligand_id": lig_id, "receptor": receptor_name, "success": False, "error": f"PDBQT: {se[:80]}"}

    # Vina
    cx, cy, cz = cfg["center"]
    sx, sy, sz = cfg["size"]
    rc, so, se = run_wsl(
        f"vina --receptor {cfg['pdbqt']} --ligand {win_to_wsl(pdbqt)} "
        f"--out {win_to_wsl(out)} --center_x {cx} --center_y {cy} --center_z {cz} "
        f"--size_x {sx} --size_y {sy} --size_z {sz} "
        f"--exhaustiveness 8 --num_modes 9 --energy_range 4", 1200)
    if rc != 0:
        return {"ligand_id": lig_id, "receptor": receptor_name, "success": False, "error": f"Vina: {se[:80]}"}

    energies = parse_vina(so)
    if not energies:
        return {"ligand_id": lig_id, "receptor": receptor_name, "success": False, "error": "No energy", "stdout": so[:200]}

    return {
        "ligand_id": lig_id, "receptor": receptor_name, "success": True,
        "best_affinity": min(energies), "all_affinities": energies, "num_modes": len(energies)
    }

def main():
    print(f"Batch Docking: 75 ligands x 2 receptors\nStarted: {datetime.now()}\n")

    # Load ligands
    ligands = list(csv.DictReader(open(LIGAND_CSV)))
    print(f"Loaded {len(ligands)} ligands")

    results = []
    total = len(ligands) * 2
    for i, row in enumerate(ligands):
        lig_id = row.get("id", f"LIG{i:04d}")
        smiles = row["smiles"]
        name = row.get("name", lig_id)

        for rec_name, cfg in RECEPTORS.items():
            r = dock_one(rec_name, cfg, lig_id, smiles, name)
            results.append(r)
            status = f"{r['best_affinity']:.1f}" if r["success"] else r["error"]
            print(f"  [{len(results)}/{total}] {lig_id} vs {rec_name}: {status}")

    # Save
    success = sum(1 for r in results if r["success"])
    print(f"\n{'='*60}")
    print(f"Complete: {success}/{len(results)} successful ({success/len(results)*100:.0f}%)")
    print(f"Time: {datetime.now()}")

    with open(OUTPUT_DIR / "docking_results.json", "w") as f:
        json.dump(results, f, indent=2)

    # Summary CSV
    with open(OUTPUT_DIR / "docking_summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ligand_id", "name", "receptor", "best_affinity_kcal_mol", "num_modes", "success"])
        for r in results:
            if r["success"]:
                w.writerow([r["ligand_id"], r.get("name",""), r["receptor"], r["best_affinity"], r["num_modes"], True])
            else:
                w.writerow([r["ligand_id"], r.get("name",""), r["receptor"], None, 0, False])

    # Top hits
    successful = [r for r in results if r["success"]]
    if successful:
        successful.sort(key=lambda r: r["best_affinity"])
        print(f"\nTop 10 hits:")
        for r in successful[:10]:
            print(f"  {r['ligand_id']} vs {r['receptor']}: {r['best_affinity']:.1f} kcal/mol")

    print(f"\nResults: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
