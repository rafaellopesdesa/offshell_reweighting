#!/usr/bin/env python3
"""Export validated resource and covariance reports without changing scan assets."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--analysis', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    resources = json.loads((args.run/'resources/public/summary.json').read_text())
    diagnostics = json.loads((args.run/'diagnostics_report/summary.json').read_text())
    family = json.loads((args.output/'data/family_snapshot.json').read_text())
    assert resources['complete'] and diagnostics['complete'] and family['complete']
    assert len(resources['models']) == 6 and set(diagnostics['samples']) == {'4096','8192'}
    assert [m['scan_csv_sha256'] for m in resources['models']] == [m['scan_csv_sha256'] for m in family['models']]
    family_by_label = {m['label']: m for m in family['models']}
    for size, sample in diagnostics['samples'].items():
        assert sample['input_sha256'] == family_by_label[f'data_like_{size}_unpruned']['workspace_sha256']
        assert len(sample['impacts']) == 130 and sum(r['constrained'] for r in sample['impacts']) == 127
        assert sample['negative_rows'] == (2 if size == '8192' else 0)
        for view in ('pulls_nominal','pulls_asimov','auxiliary_impacts','source_groups'):
            for ext in ('svg','pdf'):
                assert (args.run/'diagnostics_report'/f'data_like_{size}_{view}.{ext}').is_file()
    for model in resources['models']:
        assert model['points'] == 61 and model['computed_points'] == 66
        for ext in ('svg','pdf'):
            assert (args.run/'resources'/f"{model['label']}_resources.{ext}").is_file()
    artifacts = {}
    for directory in (args.run/'resources', args.run/'resources/public', args.run/'diagnostics_report'):
        for source in sorted(directory.iterdir()):
            if source.suffix not in ('.csv','.svg','.pdf'): continue
            if directory == args.run/'resources' and source.suffix == '.csv': continue
            target = args.output/('data' if source.suffix == '.csv' else 'assets')/('summary_'+source.name)
            shutil.copyfile(source,target)
            artifacts[str(target.relative_to(args.output))] = hashlib.sha256(target.read_bytes()).hexdigest()
    native = json.loads((args.run/'native_workspace_graph.json').read_text())
    assert native['data_rows'] == 2077 and native['channels']['SR']['count'] == 2074
    assert native['selected_nodes']['pdf_process_S']['self_normalized']
    assert native['model_config_sets']['nuisances']['count'] == 127 and native['global_observations_all_zero']
    graph_target = args.output/'data/summary_native_workspace_graph.json'
    shutil.copyfile(args.run/'native_workspace_graph.json', graph_target)
    artifacts[str(graph_target.relative_to(args.output))] = hashlib.sha256(graph_target.read_bytes()).hexdigest()
    example = args.output/'data/summary_normalizer_usage_example.py'
    shutil.copyfile(args.analysis/'normalizer_usage_example.py',example)
    artifacts[str(example.relative_to(args.output))] = hashlib.sha256(example.read_bytes()).hexdigest()
    sources = {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [args.analysis/name for name in
        ('report_summary_resources.py','report_summary_diagnostics.py','compute_observation_response.py','normalizer_usage_example.py','inspect_native_workspace.py')]}
    summary = {'complete':True,'created_utc':datetime.now(timezone.utc).isoformat(),
        'resources':resources,'diagnostics':diagnostics,'native_workspace':native,'artifact_sha256':artifacts,'analysis_source_sha256':sources,
        'historical_benchmark_scope':'One profiled mu point; 24 hours and 500 GB reported by the author. Historical CPU count, hardware and RAM measurement convention unspecified.'}
    encoded=json.dumps(summary,indent=2,allow_nan=False)+'\n'
    (args.output/'data/summary_snapshot.json').write_text(encoded)
    (args.output/'data/summary_snapshot.js').write_text('window.SUMMARY_RESULTS = '+encoded.rstrip()+';\n')
    print(f'Exported {len(artifacts)} artifacts, six resource curves and both validated covariance reports')


if __name__ == '__main__':
    main()
