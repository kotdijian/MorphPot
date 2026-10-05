import csv
import numpy as np
import pytest
from pottery_rim_interval_plot import load_comparisons, summarize_metric


def table(tmp_path,step,rows):
    path=tmp_path/f'comparison_{step:02d}deg';path.mkdir()
    with (path/'similarity_all_division_comparison.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def record(run,station,thickness='',reference='dense'):
    return dict(run=run,reference=reference,reference_profile_count=360,mode='similarity',station_mm=station,
        thickness_delta_mm=thickness,midline_distance_mm=0.2,tip_aligned_midline_distance_mm=0.1,tangent_angle_deg=1.)


def test_common_stations_and_missing_thickness_no_zero_fill(tmp_path):
    table(tmp_path,5,[record('a',5,1),record('b',5,''),record('a',10,2)])
    table(tmp_path,10,[record('c',5,3),record('c',15,9)])
    tables,stations,_=load_comparisons(tmp_path,None,'similarity','all')
    assert np.array_equal(stations,[5])
    stat=summarize_metric(tables[5],stations,'thickness_delta_mm')
    assert np.array_equal(stat,[[1,1,1,1]])


def test_signed_phase_median_and_full_range():
    rows=[record('a',5,-2),record('b',5,1),record('c',5,9)]
    assert np.array_equal(summarize_metric(rows,[5],'thickness_delta_mm'),[[1,-2,9,3]])


def test_mixed_references_rejected(tmp_path):
    table(tmp_path,5,[record('a',5,1)])
    table(tmp_path,10,[record('b',5,1,reference='other')])
    with pytest.raises(ValueError,match='same reference'):load_comparisons(tmp_path,None,'similarity','all')


def test_duplicate_phase_station_rejected(tmp_path):
    table(tmp_path,5,[record('a',5,1),record('a',5,2)])
    table(tmp_path,10,[record('b',5,1)])
    with pytest.raises(ValueError,match='Duplicate'):load_comparisons(tmp_path,None,'similarity','all')
