"""Export the existing matched-noise geometry without computing persistence.

The saved seeds, normal directions and signed draws are reused at every
amplitude. This reconstructs the inputs to the retained observations; it
does not add samples or change any released result.
"""
from pathlib import Path
import argparse
import json
import sys

sys.dont_write_bytecode = True

import numpy as np
import paper1_methods as m


def coordinate_arrays(config):
    arrays = {}
    for seed in config['seeds']:
        latent, clean = m.sample_cloud(seed, 'uniform', n=config['n'])
        # Keep the same arithmetic as the original noise trial reconstruction.
        normal = clean * (2 * np.pi)
        signed_draw = np.random.default_rng(seed + m.NOISE_SEED_OFFSET).uniform(
            -1., 1., size=config['n'])
        arrays[f'{seed}_latent'] = latent
        arrays[f'{seed}_clean'] = clean
        arrays[f'{seed}_normal'] = normal
        arrays[f'{seed}_signed_draw'] = signed_draw
        for delta in config['delta']:
            arrays[f'{seed}_delta{delta:g}_observed'] = (
                clean + delta * signed_draw[:, None] * normal)
    return arrays


def export(release, output):
    config = json.loads((release / 'results/matched_noise/config.json').read_text())
    arrays = coordinate_arrays(config)
    with output.open('wb') as stream:
        np.savez_compressed(stream, **arrays)
    print(f'Exported {len(arrays)} arrays for {len(config["seeds"])} retained noise seeds.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--release', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output', type=Path, required=True,
                        help='Scratch NPZ output outside final_release folder')
    args = parser.parse_args()
    release, output = args.release.resolve(), args.output.resolve()
    if output == release or output.is_relative_to(release):
        parser.error('Choose scratch output outside final_release')
    export(release, output)
