"""
CLI module for Protein Structure Analysis.
"""
import argparse
import json
import os
import sys
from typing import List

from .engine import (
    Residue,
    classify_ramachandran,
    analyze_ramachandran,
    ramachandran_summary,
    assign_secondary_structure,
    calculate_hydrogen_bonds,
    generate_contact_map,
    analyze_structure,
)


def _validate_input_file(filepath: str) -> None:
    """Validate that the input file exists and is readable."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Input file not found: {filepath}")
    if not os.path.isfile(filepath):
        raise ValueError(f"Path is not a file: {filepath}")
    if not os.access(filepath, os.R_OK):
        raise PermissionError(f"Cannot read file: {filepath}")


def _parse_residues_from_json(filepath: str) -> List[Residue]:
    """Load residues from a JSON file.

    Expected format:
    [
        {"name": "ALA", "index": 1, "phi": -57.0, "psi": -47.0},
        ...
    ]
    Or with coordinates:
    [
        {"name": "ALA", "index": 1, "n": [x,y,z], "ca": [x,y,z], "c": [x,y,z], "o": [x,y,z]},
        ...
    ]
    """
    _validate_input_file(filepath)

    with open(filepath, 'r') as f:
        data = json.load(f)

    if not isinstance(data, list):
        raise ValueError(f"Expected JSON array of residues, got {type(data).__name__}")

    residues = []
    for idx, item in enumerate(data):
        if not isinstance(item, dict):
            raise ValueError(f"Residue at position {idx} is not a JSON object")

        res = Residue(
            name=item.get('name', 'ALA'),
            index=item.get('index', idx),
            phi=item.get('phi'),
            psi=item.get('psi'),
            omega=item.get('omega'),
        )
        if 'n' in item:
            res.n_coord = tuple(item['n'])
        if 'ca' in item:
            res.ca_coord = tuple(item['ca'])
        if 'c' in item:
            res.c_coord = tuple(item['c'])
        if 'o' in item:
            res.o_coord = tuple(item['o'])
        residues.append(res)

    return residues


def cmd_rama(args):
    """Run Ramachandran plot analysis."""
    residues = _parse_residues_from_json(args.input)
    classifications = analyze_ramachandran(residues)
    summary = ramachandran_summary(classifications)
    
    if args.json:
        output = {
            'classifications': [
                {
                    'residue_index': c.residue_index,
                    'residue_name': c.residue_name,
                    'phi': round(c.phi, 1),
                    'psi': round(c.psi, 1),
                    'region': c.region,
                    'structure_type': c.structure_type,
                }
                for c in classifications
            ],
            'summary': summary,
        }
        print(json.dumps(output, indent=2))
    else:
        print(f"{'Res':>5} {'Name':>4} {'Phi':>8} {'Psi':>8} {'Region':>10} {'Type':>12}")
        print("-" * 55)
        for c in classifications:
            print(f"{c.residue_index:>5} {c.residue_name:>4} {c.phi:>8.1f} {c.psi:>8.1f} "
                  f"{c.region:>10} {c.structure_type:>12}")
        print(f"\nSummary:")
        print(f"  Favored:  {summary['favored']} ({summary['favored_pct']:.1f}%)")
        print(f"  Allowed:  {summary['allowed']} ({summary['allowed_pct']:.1f}%)")
        print(f"  Outlier:  {summary['outlier']} ({summary['outlier_pct']:.1f}%)")
    
    return 0


def cmd_ss(args):
    """Assign secondary structure."""
    residues = _parse_residues_from_json(args.input)
    assignments = assign_secondary_structure(residues)
    
    if args.json:
        output = [
            {'residue_index': residues[i].index, 'residue_name': residues[i].name, 'ss': ss}
            for i, ss in enumerate(assignments)
        ]
        print(json.dumps(output, indent=2))
    else:
        ss_names = {'H': 'Helix', 'E': 'Strand', 'T': 'Turn', 'P': 'PPII', 'C': 'Coil'}
        print(f"{'Res':>5} {'Name':>4} {'SS':>3} {'Type':>8}")
        print("-" * 25)
        for i, ss in enumerate(assignments):
            print(f"{residues[i].index:>5} {residues[i].name:>4} {ss:>3} {ss_names.get(ss, 'Unknown'):>8}")
        
        counts = {}
        for ss in assignments:
            counts[ss] = counts.get(ss, 0) + 1
        print(f"\nDistribution:")
        for ss, count in sorted(counts.items()):
            print(f"  {ss_names.get(ss, ss)}: {count}")
    
    return 0


def cmd_hbonds(args):
    """Calculate hydrogen bonds."""
    residues = _parse_residues_from_json(args.input)
    hbonds = calculate_hydrogen_bonds(residues, args.cutoff)
    
    if args.json:
        output = [
            {
                'donor': h.donor_residue,
                'acceptor': h.acceptor_residue,
                'distance': h.distance,
                'energy': h.energy,
            }
            for h in hbonds
        ]
        print(json.dumps(output, indent=2))
    else:
        print(f"Found {len(hbonds)} hydrogen bonds (cutoff: {args.cutoff} A)")
        print(f"{'Donor':>6} {'Acceptor':>8} {'Dist(A)':>8} {'E(kcal)':>8}")
        print("-" * 35)
        for h in hbonds:
            print(f"{h.donor_residue:>6} {h.acceptor_residue:>8} {h.distance:>8.3f} {h.energy:>8.3f}")
    
    return 0


def cmd_contacts(args):
    """Generate contact map."""
    residues = _parse_residues_from_json(args.input)
    cmap = generate_contact_map(residues, args.cutoff)
    
    if args.json:
        output = {
            'size': cmap.size,
            'cutoff': cmap.cutoff,
            'num_contacts': len(cmap.contacts),
            'contacts': [{'residue_i': c[0], 'residue_j': c[1]} for c in cmap.contacts],
        }
        print(json.dumps(output, indent=2))
    else:
        print(f"Contact map: {cmap.size} residues, cutoff {cmap.cutoff} A")
        print(f"Total contacts: {len(cmap.contacts)}")
        for c in cmap.contacts[:20]:
            print(f"  {c[0]} - {c[1]}")
        if len(cmap.contacts) > 20:
            print(f"  ... and {len(cmap.contacts) - 20} more")
    
    return 0


def cmd_analyze(args):
    """Run full structure analysis."""
    residues = _parse_residues_from_json(args.input)
    result = analyze_structure(residues, args.hbond_cutoff, args.contact_cutoff)
    
    if args.json:
        output = {
            'num_residues': result.num_residues,
            'secondary_structure': result.secondary_structure,
            'ramachandran_summary': result.summary.get('ramachandran', {}),
            'num_hbonds': result.num_hbonds if 'num_hbonds' in result.summary else len(result.hbonds),
            'num_contacts': result.num_contacts if 'num_contacts' in result.summary else (
                len(result.contact_map.contacts) if result.contact_map else 0
            ),
            'summary': result.summary,
        }
        print(json.dumps(output, indent=2))
    else:
        print(f"Protein Structure Analysis")
        print(f"{'=' * 50}")
        print(f"Residues: {result.num_residues}")
        print(f"\nSecondary Structure:")
        ss_names = {'H': 'Helix', 'E': 'Strand', 'T': 'Turn', 'P': 'PPII', 'C': 'Coil'}
        for ss, count in result.summary.get('secondary_structure_counts', {}).items():
            print(f"  {ss_names.get(ss, ss)}: {count}")
        print(f"\nRamachandran:")
        rama = result.summary.get('ramachandran', {})
        print(f"  Favored: {rama.get('favored', 0)} ({rama.get('favored_pct', 0):.1f}%)")
        print(f"  Allowed: {rama.get('allowed', 0)} ({rama.get('allowed_pct', 0):.1f}%)")
        print(f"  Outlier: {rama.get('outlier', 0)} ({rama.get('outlier_pct', 0):.1f}%)")
        print(f"\nHydrogen bonds: {len(result.hbonds)}")
        if result.contact_map:
            print(f"Contacts: {len(result.contact_map.contacts)}")
    
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog='protein-folding-alpha-geometry',
        description='Protein Structure Analysis: Ramachandran plots, secondary structure, H-bonds, contact maps',
    )
    subparsers = parser.add_subparsers(dest='command', required=True)
    
    # Ramachandran
    p_rama = subparsers.add_parser('rama', help='Ramachandran plot analysis')
    p_rama.add_argument('-i', '--input', required=True, help='JSON file with residue angles')
    p_rama.add_argument('--json', action='store_true', help='Output as JSON')
    
    # Secondary structure
    p_ss = subparsers.add_parser('ss', help='Secondary structure assignment')
    p_ss.add_argument('-i', '--input', required=True, help='JSON file with residue angles')
    p_ss.add_argument('--json', action='store_true', help='Output as JSON')
    
    # Hydrogen bonds
    p_hbonds = subparsers.add_parser('hbonds', help='Hydrogen bond calculation')
    p_hbonds.add_argument('-i', '--input', required=True, help='JSON file with coordinates')
    p_hbonds.add_argument('--cutoff', type=float, default=3.5, help='N-O distance cutoff (A)')
    p_hbonds.add_argument('--json', action='store_true', help='Output as JSON')
    
    # Contact map
    p_contacts = subparsers.add_parser('contacts', help='Contact map generation')
    p_contacts.add_argument('-i', '--input', required=True, help='JSON file with coordinates')
    p_contacts.add_argument('--cutoff', type=float, default=8.0, help='CA-CA distance cutoff (A)')
    p_contacts.add_argument('--json', action='store_true', help='Output as JSON')
    
    # Full analysis
    p_analyze = subparsers.add_parser('analyze', help='Full structure analysis')
    p_analyze.add_argument('-i', '--input', required=True, help='JSON file with residue data')
    p_analyze.add_argument('--hbond-cutoff', type=float, default=3.5, help='H-bond cutoff (A)')
    p_analyze.add_argument('--contact-cutoff', type=float, default=8.0, help='Contact cutoff (A)')
    p_analyze.add_argument('--json', action='store_true', help='Output as JSON')
    
    args = parser.parse_args(argv)

    commands = {
        'rama': cmd_rama,
        'ss': cmd_ss,
        'hbonds': cmd_hbonds,
        'contacts': cmd_contacts,
        'analyze': cmd_analyze,
    }

    try:
        return commands[args.command](args)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON in input file - {e}", file=sys.stderr)
        return 1
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    except PermissionError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
