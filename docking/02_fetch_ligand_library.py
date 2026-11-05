import os
import time
import logging
from typing import List, Dict, Tuple

import requests
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from rdkit import DataStructs
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

OUT_DIR = os.path.join('results', 'molecular_docking')
LIGAND_CSV = os.path.join('ligands', '01_candidate_ligand_library.csv')
SUMMARY_CSV = os.path.join(OUT_DIR, '01_ligand_library_summary.csv')

CHEMBL_BASE = 'https://www.ebi.ac.uk/chembl/api/data/molecule'
HEADERS = {'User-Agent': 'HPV-E6-Library-Build/1.1'}

QUERIES = [
    'antiviral', 'antineoplastic', 'oncology', 'antiinfective', 'HPV', 'papillomavirus',
    'cancer', 'tumor', 'neoplasm', 'chemotherapy', 'oncolytic', 'antitumor', 'antineoplastic agents',
    'virus inhibitor', 'viral', 'DNA virus', 'E6', 'p53', 'PDZ'
]

def ensure_dirs():
    os.makedirs('ligands', exist_ok=True)
    os.makedirs(OUT_DIR, exist_ok=True)

def fetch_chembl_search(query: str, page_size: int = 200, max_pages: int = 100) -> List[Dict]:
    results = []
    for page in range(0, max_pages):
        params = {'q': query, 'format': 'json', 'page_size': page_size, 'page': page+1}
        try:
            r = requests.get(f'{CHEMBL_BASE}/search', params=params, headers=HEADERS, timeout=30)
            if r.status_code != 200:
                logging.warning(f'ChEMBL search failed: {query} page {page+1} status {r.status_code}')
                break
            data = r.json()
            items = data.get('molecules', [])
            if not items:
                break
            for it in items:
                mstruct = it.get('molecule_structures') or {}
                smiles = mstruct.get('canonical_smiles')
                chembl_id = it.get('molecule_chembl_id')
                name = it.get('pref_name') or chembl_id
                if smiles:
                    results.append({'source': 'ChEMBL', 'id': chembl_id, 'name': name, 'smiles': smiles})
            if len(items) < page_size:
                break
            time.sleep(0.1)
        except Exception as e:
            logging.error(f'Error fetching ChEMBL for query {query}: {e}')
            break
    return results

def compute_props(smiles: str) -> Dict:
    m = Chem.MolFromSmiles(smiles)
    if not m:
        return {}
    return {
        'MW': Descriptors.MolWt(m),
        'LogP': Crippen.MolLogP(m),
        'HBD': rdMolDescriptors.CalcNumHBD(m),
        'HBA': rdMolDescriptors.CalcNumHBA(m),
        'RotBonds': rdMolDescriptors.CalcNumRotatableBonds(m),
        'TPSA': rdMolDescriptors.CalcTPSA(m),
        'Rings': rdMolDescriptors.CalcNumRings(m)
    }

def lipinski_veber_filter(props: Dict) -> Tuple[bool, bool]:
    if not props:
        return False, False
    lipinski = (props['MW'] <= 500 and props['LogP'] <= 5 and props['HBD'] <= 5 and props['HBA'] <= 10)
    veber = (props['RotBonds'] <= 10 and props['TPSA'] <= 140)
    return lipinski, veber

def smiles_fingerprint(smiles: str):
    m = Chem.MolFromSmiles(smiles)
    return rdMolDescriptors.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) if m else None

def diversity_filter(rows: List[Dict], threshold: float = 0.8, target_count: int = 1000) -> List[Dict]:
    selected = []
    fps = []
    for r in rows:
        fp = smiles_fingerprint(r['smiles'])
        if fp is None:
            continue
        if not fps:
            selected.append(r)
            fps.append(fp)
        else:
            sims = [DataStructs.TanimotoSimilarity(fp, xfp) for xfp in fps]
            if max(sims) < threshold:
                selected.append(r)
                fps.append(fp)
        if len(selected) >= target_count:
            break
    return selected

def build_library(target_count: int = 1000) -> pd.DataFrame:
    ensure_dirs()
    all_rows: List[Dict] = []
    for q in QUERIES:
        logging.info(f'Fetching from ChEMBL: query={q}')
        rows = fetch_chembl_search(q, page_size=200, max_pages=100)
        all_rows.extend(rows)
        logging.info(f'ChEMBL fetched {len(rows)} for query {q}, total {len(all_rows)}')
    df = pd.DataFrame(all_rows)
    df = df.dropna(subset=['smiles']).drop_duplicates(subset=['smiles'])
    logging.info(f'Unique SMILES count after de-dup: {len(df)}')
    props_list = []
    for smi in tqdm(df['smiles'].tolist(), desc='Computing properties'):
        props_list.append(compute_props(smi))
    props_df = pd.DataFrame(props_list)
    df = pd.concat([df.reset_index(drop=True), props_df], axis=1)
    df['LipinskiPass'], df['VeberPass'] = zip(*df.apply(lambda r: lipinski_veber_filter({
        'MW': r['MW'], 'LogP': r['LogP'], 'HBD': r['HBD'], 'HBA': r['HBA'], 'RotBonds': r['RotBonds'], 'TPSA': r['TPSA']
    }), axis=1))
    filtered = df[(df['LipinskiPass']) & (df['VeberPass'])]
    logging.info(f'After Lipinski+Veber filter: {len(filtered)} molecules')
    rows_dicts = filtered.to_dict(orient='records')
    diverse = diversity_filter(rows_dicts, threshold=0.8, target_count=target_count)
    diverse_df = pd.DataFrame(diverse)
    if len(diverse_df) < target_count:
        more = diversity_filter(rows_dicts, threshold=0.9, target_count=target_count)
        diverse_df = pd.concat([diverse_df, pd.DataFrame(more)]).drop_duplicates(subset=['smiles'])
    logging.info(f'Final candidate count: {len(diverse_df)}')
    diverse_df = diverse_df[['source','id','name','smiles','MW','LogP','HBD','HBA','RotBonds','TPSA','Rings','LipinskiPass','VeberPass']]
    os.makedirs(os.path.dirname(LIGAND_CSV), exist_ok=True)
    diverse_df.to_csv(LIGAND_CSV, index=False)
    summary = {
        'total_initial': int(len(df)),
        'after_filters': int(len(filtered)),
        'final_candidates': int(len(diverse_df))
    }
    pd.DataFrame([summary]).to_csv(SUMMARY_CSV, index=False)
    logging.info(f'Wrote {LIGAND_CSV} and {SUMMARY_CSV}')
    return diverse_df

if __name__ == '__main__':
    build_library(target_count=1000)