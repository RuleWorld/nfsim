"""Focused XML regressions for RuleWorld/nfsim#102, #104, #106, #107 and #115."""
import pathlib
import math
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

BINARY = str(pathlib.Path(sys.argv[1]).resolve())
FIXTURES = pathlib.Path(__file__).resolve().parent


def run(root, folder, name, extra=(), seed=123):
    xml = folder / (name + '.xml')
    ET.ElementTree(root).write(xml)
    output = folder / (name + '.gdat')
    proc = subprocess.run([BINARY, '-xml', str(xml), '-o', str(output),
                           '-sim', '2', '-oSteps', '4', '-seed', str(seed), *extra],
                          capture_output=True, text=True, timeout=30)
    text = proc.stdout + proc.stderr
    return proc, text, output.read_text() if output.exists() else ''


def final_mean(root, folder, name, column, seeds=200):
    values = []
    for seed in range(seeds):
        _, log, data = run(root, folder, '{0}_{1}'.format(name, seed),
                           seed=123 + seed)
        rows = [line.split() for line in data.splitlines()
                if line and not line.startswith('#')]
        header = next(line.lstrip('# ').split() for line in data.splitlines()
                      if line.startswith('#'))
        assert rows and column in header, (name, column, data, log)
        values.append(float(rows[-1][header.index(column)]))
    return sum(values) / len(values)


def load(name):
    root = ET.parse(FIXTURES / (name + '.xml')).getroot()
    for node in root.iter():
        node.tag = node.tag.split('}')[-1]
    return root


with tempfile.TemporaryDirectory() as temporary:
    folder = pathlib.Path(temporary)
    mixed = load('issue104')
    mixed_ab = final_mean(mixed, folder, 'mixed', 'MAB')
    # A local function that also reads a system observable outside its
    # argument must follow that observable. The published BNG reference for
    # this model at t=2 is AB ~= 29.8; freezing MAB at its initial value (the
    # old stale behaviour) gives ~26.4.
    assert 28.5 < mixed_ab < 31.0, mixed_ab

    local = load('issue104')
    # Keep the scoped term while removing the outside observable dependency.
    function = local.find('.//Function')
    expr = function.find('Expression')
    expr.text = 'k*(1+AC(x))'
    refs = function.find('ListOfReferences')
    for ref in list(refs):
        if ref.attrib.get('name') == 'MAB':
            refs.remove(ref)
    local_ab = final_mean(local, folder, 'local', 'MAB')
    assert local_ab > 24.0, local_ab
    # The stale value is exactly this no-global-term control, so the mixed
    # model must sit clearly above it.
    assert mixed_ab > local_ab + 2.0, (mixed_ab, local_ab)

    pattern = load('issue106')
    _, log, data = run(pattern, folder, 'pattern')
    assert data, log
    rows = [line.split() for line in data.splitlines() if line and not line.startswith('#')]
    assert list(map(float, rows[0][1:])) == [30, 20, 20, 10], rows[0]
    assert all(float(row[4]) == 10 for row in rows), rows
    molecules = load('issue106')
    for parent in molecules.iter():
        if parent.tag in ('Pattern', 'ReactantPattern', 'ProductPattern'):
            compartment = parent.attrib.pop('compartment', None)
            if compartment:
                for molecule in parent.findall('./ListOfMolecules/Molecule'):
                    molecule.attrib.setdefault('compartment', compartment)
    _, log, control = run(molecules, folder, 'molecules')
    assert data == control, (data, control, log)
    _, log, connected = run(pattern, folder, 'connected', ('-connect',))
    assert data == connected, (data, connected, log)

    override = load('issue106')
    observable = override.find(".//Observable[@name='AECp']/ListOfPatterns/Pattern")
    observable.find('./ListOfMolecules/Molecule').set('compartment', 'CP')
    _, log, data_override = run(override, folder, 'override')
    first = [line.split() for line in data_override.splitlines() if line and not line.startswith('#')][0]
    assert float(first[3]) == 10, first

    sparse = load('issue107')
    _, log, off = run(sparse, folder, 'sparse_off')
    _, log, on = run(sparse, folder, 'sparse_on', ('-connect',))
    assert off and off == on, (off, on, log)

    for tag in ('Pattern', 'ReactantPattern', 'ProductPattern'):
        invalid = load('issue106')
        invalid.find('.//' + tag).set('compartment', 'missing')
        _, log, _ = run(invalid, folder, 'unknown_' + tag)
        assert "unknown compartment 'missing'" in log, log

    # Numeric Arrhenius rates, directionality, expression-valued energies,
    # and RT=1 energy units (RT remains an ordinary XML parameter).
    energy = load('issue102')

    def expanded_rates(model, name):
        xml = folder / (name + '.xml')
        ET.ElementTree(model).write(xml)
        proc = subprocess.run([BINARY, '-xml', str(xml), '-o', str(folder / (name + '.gdat')),
                               '-v', '-sim', '0.001',
                               '-oSteps', '1', '-seed', '1'],
                              capture_output=True, text=True, timeout=30)
        log = proc.stdout + proc.stderr
        assert proc.returncode == 0, log
        rates = {}
        for direction, rule, value in re.findall(
                r'Created (forward|reverse)(?: state-change)? rule: (\S+)\s+rate=([0-9.eE+-]+)', log):
            rates[rule] = float(value)
        return rates

    rates_25 = expanded_rates(energy, 'energy_rt25')
    expected = {
        'Reversible_fwd': math.exp(-0.25 * math.log(1000)),
        '_reverse_Reversible_rev': math.exp(0.75 * math.log(1000)),
        'Unbind_rev': math.exp(-0.75 * math.log(1000)),
        'Bind_fwd': 1e-5,
        'Flip_fwd': math.exp(-0.25 * math.log(1000) - math.log(1000)),
        '_reverse_Flip_fwd': math.exp(0.75 * math.log(1000) - math.log(1000)),
        'OneWayUp_fwd': math.exp(-0.25 * math.log(1000)),
        'OneWayDown_fwd': math.exp(0.25 * math.log(1000)),
        'IndependentBind_fwd': math.exp(-0.25 * math.log(1000)),
        'IndependentUnbind_rev': math.exp(0.25 * math.log(1000)),
    }
    assert rates_25.keys() == expected.keys(), rates_25
    for rule, rate in expected.items():
        assert math.isclose(rates_25[rule], rate, rel_tol=1e-5), (rule, rates_25[rule], rate)

    # RT and dG are both rescaled, so expression values and all rates remain
    # identical; no extra RT division is applied by NFsim.
    energy_rt1 = load('issue102')
    for parameter in energy_rt1.findall('.//ListOfParameters/Parameter'):
        if parameter.attrib.get('id') == 'RT':
            parameter.set('value', '1')
            parameter.set('expr', '1')
        elif parameter.attrib.get('id') == 'dG':
            parameter.set('value', str(math.log(1000)))
            parameter.set('expr', str(math.log(1000)))
    rates_1 = expanded_rates(energy_rt1, 'energy_rt1')
    assert rates_1.keys() == rates_25.keys(), (rates_1, rates_25)
    for rule in rates_25:
        assert math.isclose(rates_1[rule], rates_25[rule], rel_tol=1e-5), (
            rule, rates_1[rule], rates_25[rule])

    for expression, expected_error in (('missing_energy', 'cannot evaluate energy expression'),
                                       ('1/0', 'not finite')):
        invalid = load('issue102')
        invalid.find('.//EnergyPattern').set('expression', expression)
        path = folder / ('invalid_energy_' + expression.replace('/', '_') + '.xml')
        ET.ElementTree(invalid).write(path)
        proc = subprocess.run([BINARY, '-xml', str(path), '-o', str(folder / (path.stem + '.gdat')), '-v'],
                              capture_output=True, text=True, timeout=30)
        assert expected_error in proc.stdout + proc.stderr, proc.stdout + proc.stderr

    overflow = load('issue102')
    overflow.find(".//Parameter[@id='Ea0']").set('value', '-1000')
    overflow.find(".//Parameter[@id='Ea0']").set('expr', '-1000')
    overflow_path = folder / 'overflow_rate.xml'
    ET.ElementTree(overflow).write(overflow_path)
    proc = subprocess.run([BINARY, '-xml', str(overflow_path), '-o', str(folder / 'overflow_rate.gdat'), '-v'],
                          capture_output=True, text=True, timeout=30)
    assert 'not finite and non-negative' in proc.stdout + proc.stderr, proc.stdout + proc.stderr

    # Native NFsim parity cases for issue #109 and #110, generated by BNG2.
    inheritance = load('issue109')
    _, log, data = run(inheritance, folder, 'issue109')
    assert data, log
    rows = [line.split() for line in data.splitlines() if line and not line.startswith('#')]
    assert all(float(rows[-1][index]) > 0 for index in (1, 2, 3)), rows[-1]

    placement = load('issue110_seed_bind')
    _, log, data = run(placement, folder, 'issue110_seed_bind')
    assert data, log
    rows = [line.split() for line in data.splitlines() if line and not line.startswith('#')]
    assert [float(rows[0][index]) for index in (1, 2)] == [100, 0], rows[0]
    assert float(rows[-1][3]) > 0, rows[-1]
    assert all(float(row[4]) == 0 for row in rows), rows

    connected_move = load('issue110_move')
    _, log, data = run(connected_move, folder, 'issue110_move')
    assert data, log
    rows = [line.split() for line in data.splitlines() if line and not line.startswith('#')]
    assert float(rows[-1][1]) == 100 and float(rows[-1][2]) == 0, rows[-1]
    assert float(rows[-1][3]) > 0 and float(rows[-1][4]) > 0, rows[-1]

    # Issue #115. A firing takes its product molecules out of their observables
    # and puts them back; a rate law read in between is read at a state the
    # system is never in. With one X, kgap*(1+9*Xp/Xtot) was 0/0 there: a_tot
    # became NaN and the run ended silently at the first X firing (the last
    # change in 20 seeds came at t = 1 to 62 of 500). The second model deletes
    # a molecule in the same firing, which updates propensities by another path
    # than the observables. In the third, 1+9*(Xtot-Xp) is negative in between
    # and in no state the system reaches, which ended the run with a crash. In
    # the fourth, kd/Xp with TotalRate is infinite once Xp is counted again and
    # before X has left the rule's reactants; in the fifth the same rule has two
    # reactants, and the rate was read with one of them back in the rule and
    # the other not.
    def issue115_runs(name, seeds, extra=()):
        xml = folder / (name + '.xml')
        ET.ElementTree(load(name)).write(xml)
        for seed in seeds:
            output = folder / '{0}_{1}.gdat'.format(name, seed)
            proc = subprocess.run([BINARY, '-xml', str(xml), '-o', str(output), '-sim', '500',
                                   '-oSteps', '500', '-seed', str(seed), *extra],
                                  capture_output=True, text=True, timeout=120)
            assert proc.returncode == 0 and output.exists(), (name, seed, proc.stdout + proc.stderr)
            yield seed, [[float(field) for field in line.split()]
                         for line in output.read_text().splitlines()
                         if line and not line.startswith('#')]

    # With -connect the products are put back in their reactions by another
    # function, which reads the rates too.
    for name in ('issue115', 'issue115_delete', 'issue115_negative', 'issue115_total_rate',
                 'issue115_bimolecular'):
        for extra in ((), ('-connect',)):
            for seed, rows in issue115_runs(name, range(1, 11), extra):
                changes = [after[0] for before, after in zip(rows, rows[1:])
                           if before[1:] != after[1:]]
                assert changes and changes[-1] >= 400, (
                    name, extra, seed, changes[-1] if changes else None)

    # What is held during a firing is applied after it. Here a synthesis rate,
    # kmake*Xp/Xtot, changes only when Xp does, and nothing else recomputes it:
    # if the held update were dropped, the selector would keep a stale total and
    # X would sit in its P state three quarters of the time (measured, 0.75 to
    # 0.77 over three sets of 20 seeds), where kp = kd puts it there half the
    # time (0.46 to 0.51).
    rows = [row for _, run in issue115_runs('issue115_refresh', range(1, 21)) for row in run]
    in_p = sum(row[1] for row in rows) / len(rows)
    assert 0.38 < in_p < 0.62, in_p

    # A rate law negative in a state the model reaches is refused, and the
    # message shows the offending function. For a composite rate law
    # (kgap*boost()) that function is a CompositeFunction and the plain one is
    # null: printing the plain one crashed the process (exit 139) before the
    # refusal could be made.
    negative = load('negative_composite_rate')
    xml = folder / 'negative_composite_rate.xml'
    ET.ElementTree(negative).write(xml)
    proc = subprocess.run([BINARY, '-xml', str(xml), '-o', str(folder / 'negative_composite_rate.gdat'),
                           '-sim', '10', '-oSteps', '10', '-seed', '1'],
                          capture_output=True, text=True, timeout=60)
    text = proc.stdout + proc.stderr
    assert proc.returncode == 1, (proc.returncode, text[-2000:])
    assert 'negative propensity (here -8)' in text and "Composite Function: '_rateLaw1()'" in text, text[-2000:]
    assert 'Observable:  Xtot = 1' in text, text[-2000:]

print('XML regressions passed: Arrhenius numeric rates/directions, RT rescaling, invalid/nonfinite energy and rate refusal, '
      'mixed-scope refresh, compartment inheritance/override, -connect parity, unknown compartments, and a rate law '
      'dividing by an observable a firing empties, and a negative composite rate law refused without a crash.')
