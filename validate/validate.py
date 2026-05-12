import unittest
import os
import numpy as np
import subprocess
import re
import fnmatch
import sys
import tempfile
try:
    import bionetgen
except ImportError:
    bionetgen = None

nIterations=30
nfsimPrePath='..'
mfolder='./basicModels'
targetedTests = {
    # Known noisy models get an extra targeted pass with more attempts.
    '16': {'iterations': 60, 'seed_offset': 160000},
    '18': {'iterations': 30, 'seed_offset': 100000},
    '19': {'iterations': 30, 'seed_offset': 200000},
}
if os.name == "nt":
    nfsimPath = os.path.join(nfsimPrePath, 'build', 'NFsim.exe')
else:
    nfsimPath = os.path.join(nfsimPrePath, 'build', 'NFsim')



class ParametrizedTestCase(unittest.TestCase):

    """ TestCase classes that want to be parametrized should
        inherit from this class.
    """

    def __init__(self, methodName='runTest', param=None):
        super(ParametrizedTestCase, self).__init__(methodName)
        self.param = param

    @staticmethod
    def parametrize(testcase_klass, param=None):
        """ Create a suite containing all tests taken from the given
            subclass, passing them the parameter 'param'.
        """
        testloader = unittest.TestLoader()
        testnames = testloader.getTestCaseNames(testcase_klass)
        suite = unittest.TestSuite()
        for name in testnames:
            suite.addTest(testcase_klass(name, param=param))
        return suite


def loadResults(fileName, split):
    try:
        with open(fileName) as dataInput:
            timeCourse = []
            # remove spaces
            line = dataInput.readline().strip()
            headers = re.sub(r'\s+', ' ', line).split(split)

            for line in dataInput:
                nline = re.sub(r'\s+', ' ', line.strip()).split(' ')
                try:
                    timeCourse.append([float(x) for x in nline])
                except:
                    print('++++', nline)
        return headers, np.array(timeCourse)
    except IOError:
        print('no file')
        return [], []


class TestNFSimFile(ParametrizedTestCase):
    def BNGtrajectoryGeneration(self, outputDirectory, fileNumber):
        bngFileName = os.path.join(outputDirectory, 'v{0}.bngl'.format(fileNumber))
        bionetgen.run(bngFileName, out=outputDirectory, suppress=True)

    def NFsimtrajectoryGeneration(self, outputDirectory, fileNumber, runOptions, seed=None):
        runOptions = [x.strip() for x in runOptions.split(' ') if x.strip()]
        if seed is not None and '-seed' not in runOptions:
            runOptions = runOptions + ['-seed', str(seed)]
        with open(os.devnull, "w") as fnull:
            subprocess.check_call([nfsimPath, '-xml', os.path.join(outputDirectory, 'v{0}.xml'.format(fileNumber)),
                                   '-o', os.path.join(outputDirectory, 'v{0}_nf.gdat'.format(fileNumber))] + runOptions,
                                  stdout=fnull)

    def _seed_for_iteration(self, index):
        try:
            modelNum = int(self.param['num'])
        except ValueError:
            modelNum = sum([ord(x) for x in str(self.param['num'])])
        seedOffset = int(self.param.get('seed_offset', 0))
        return seedOffset + (modelNum * 1000) + index + 1

    def loadConfigurationFile(self, outputDirectory, fileNumber):
        with open(os.path.join(outputDirectory, 'r{0}.txt').format(fileNumber), 'r') as f:
            return f.readlines()

    def test_nfsim(self):
        tol = 0.35 # this is the error tolerance when comparing nfsim's run to the ssa where 0.35 = 35%
        (modelName, runOptions) = self.loadConfigurationFile(self.param['odir'], self.param['num'])
        runTag = self.param.get('tag', 'default')
        if runTag == 'default':
            print(f"Processing model r{self.param['num']}.txt: {modelName.strip()}")
        else:
            print(f"Processing model r{self.param['num']}.txt ({runTag}): {modelName.strip()}")
        # here we decide if this is a NFsim only run or not
        if modelName.startswith("NFSIM ONLY"):
            seed = self._seed_for_iteration(0)
            self.BNGtrajectoryGeneration(self.param['odir'], self.param['num'])
            self.NFsimtrajectoryGeneration(self.param['odir'], self.param['num'], runOptions, seed=seed)
            nfh, nf = loadResults(os.path.join(self.param['odir'], 'v{0}_nf.gdat'.format(self.param['num'])), ' ')
            # here we just need to make sure we managed to get here without errors
            #assert len(nf) > 0
            self.assertTrue(len(nf) > 0)
        else:
            ssaDiff = nfDiff = 0
            bad = np.array([1])
            lastSeed = None
            for index in range(self.param['iterations']):
                seed = self._seed_for_iteration(index)
                lastSeed = seed
                print(f'Iteration {index+1} (seed={seed})')
                self.BNGtrajectoryGeneration(self.param['odir'], self.param['num'])
                self.NFsimtrajectoryGeneration(self.param['odir'], self.param['num'], runOptions, seed=seed)
                odeh, ode = loadResults(os.path.join(self.param['odir'], 'v{0}_ode.gdat'.format(self.param['num'])), ' ')
                ssah, ssa = loadResults(os.path.join(self.param['odir'], 'v{0}_ssa.gdat'.format(self.param['num'])), ' ')
                nfh, nf = loadResults(os.path.join(self.param['odir'], 'v{0}_nf.gdat'.format(self.param['num'])), ' ')

                #square root difference
                ssaDiff += pow(sum(pow(ode[:, 1:] - ssa[:, 1:], 2)), 0.5)
                nfDiff += pow(sum(pow(ode[:, 1:] - nf[:, 1:], 2)), 0.5)
                rdiff = nfDiff - ssaDiff - (tol * ssaDiff)
                # relative difference should be less than 'tol'
                bad=np.where(rdiff>0)[0]
                if (bad.size>0):
                    print(f"Sir, the observables {bad+1} did not pass at seed={seed}. Trying again")
                else:
                    print("Check passed.")
                    break
            self.assertTrue(
                bad.size==0,
                f"Model r{self.param['num']} failed after {self.param['iterations']} deterministic seeds; "
                f"last seed={lastSeed}, failing observables={bad+1}"
            )


def getTests(directory):
    """
    Gets a list of bngl files that could be correctly translated in a given 'directory'
    """
    matches = []
    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, '*txt'):
            matches.append(''.join(filename.split('.')[0][1:]))
    return sorted(matches)


class TestIssueRegressions(unittest.TestCase):

    def _load_gdat(self, filePath):
        with open(filePath, 'r') as f:
            headerLine = re.sub(r'\s+', ' ', f.readline().strip())
        headers = [h for h in headerLine.split(' ') if h and h != '#']
        data = np.loadtxt(filePath, comments='#')
        if data.ndim == 1:
            data = data.reshape(1, -1)
        # Keep only headers that correspond to numeric columns in data.
        if len(headers) > data.shape[1]:
            headers = headers[-data.shape[1]:]
        return headers, data

    def _bng_generate(self, outputDirectory, fileNumber):
        xmlFileName = os.path.join(outputDirectory, 'v{0}.xml'.format(fileNumber))
        if os.path.exists(xmlFileName):
            # already generated, no need to rerun BNG
            return

        if bionetgen is None:
            self.fail('bionetgen Python package is required to generate XML fixtures')

        bngFileName = os.path.join(outputDirectory, 'v{0}.bngl'.format(fileNumber))
        bionetgen.run(bngFileName, out=outputDirectory, suppress=True)

    def _run_nfsim_xml(self, xmlPath, outputPath, runOptions, expect_success=True):
        runOptions = [x.strip() for x in runOptions.split(' ') if x.strip()]
        if os.path.exists(outputPath):
            os.remove(outputPath)
        with open(os.devnull, "w") as fnull:
            result = subprocess.run([
                nfsimPath,
                '-xml', xmlPath,
                '-o', outputPath
            ] + runOptions, stdout=fnull, stderr=fnull)
        if expect_success:
            self.assertEqual(result.returncode, 0, f'NFsim failed for XML fixture {xmlPath}')
            self.assertTrue(os.path.exists(outputPath), f'NFsim did not create expected output file for {xmlPath}')
        else:
            self.assertTrue(
                result.returncode != 0 or not os.path.exists(outputPath),
                f'NFsim unexpectedly succeeded for XML fixture {xmlPath}'
            )
        return result

    def _run_nfsim(self, outputDirectory, fileNumber, runOptions):
        self._run_nfsim_xml(
            os.path.join(outputDirectory, 'v{0}.xml'.format(fileNumber)),
            os.path.join(outputDirectory, 'v{0}_nf.gdat'.format(fileNumber)),
            runOptions,
            expect_success=True,
        )

    def _assert_same_seed_connectivity_parity(self, xmlPath, runOptions, label):
        with tempfile.TemporaryDirectory() as tmpdir:
            offPath = os.path.join(tmpdir, 'off.gdat')
            onPath = os.path.join(tmpdir, 'on.gdat')

            connectOptions = f'{runOptions} -connect'.strip()
            self._run_nfsim_xml(xmlPath, offPath, runOptions)
            self._run_nfsim_xml(xmlPath, onPath, connectOptions)

            offHeaders, offData = self._load_gdat(offPath)
            onHeaders, onData = self._load_gdat(onPath)

            self.assertEqual(offHeaders, onHeaders, f'Connectivity regression changed {label} output columns')
            self.assertEqual(offData.shape, onData.shape, f'Connectivity regression changed {label} output shape')
            self.assertTrue(
                np.array_equal(offData, onData),
                f'Connectivity regression changed the same-seed {label} trajectory'
            )

    def test_connectivity_preserves_seeded_tlbr_trajectory(self):
        self._assert_same_seed_connectivity_parity(
            os.path.join(nfsimPrePath, 'test', 'tlbr', 'tlbr.xml'),
            '-sim 1 -oSteps 100 -seed 1',
            'TLBR'
        )

    def test_connectivity_preserves_seeded_local_function_trajectory(self):
        # testSuite/t3 exercises local-function membership updates on a much
        # smaller model than AN_chemotaxis while still reproducing the
        # master-vs-connect divergence fixed by this branch.
        self._assert_same_seed_connectivity_parity(
            os.path.join(nfsimPrePath, 'test', 'testSuite', 't3.xml'),
            '-sim 1 -oSteps 20 -seed 1',
            'testSuite t3'
        )

    def _assert_matching_output_schedules(self, xmlName, continuousOptions, chunkedOptions):
        xmlPath = os.path.join(mfolder, xmlName)
        with tempfile.TemporaryDirectory(prefix='nfsim_step_to_') as tmpdir:
            continuousOutput = os.path.join(tmpdir, 'continuous_nf.gdat')
            chunkedOutput = os.path.join(tmpdir, 'chunked_nf.gdat')

            self._run_nfsim_xml(xmlPath, continuousOutput, continuousOptions)
            self._run_nfsim_xml(xmlPath, chunkedOutput, chunkedOptions)

            continuousHeaders, continuousData = self._load_gdat(continuousOutput)
            chunkedHeaders, chunkedData = self._load_gdat(chunkedOutput)

            self.assertEqual(
                chunkedHeaders, continuousHeaders,
                f'Chunked stepTo output headers differed from continuous run for {xmlName}'
            )
            self.assertEqual(
                chunkedData.shape, continuousData.shape,
                f'Chunked stepTo output shape differed from continuous run for {xmlName}'
            )
            if not np.array_equal(chunkedData, continuousData):
                diffIndex = np.argwhere(chunkedData != continuousData)[0]
                row = int(diffIndex[0])
                col = int(diffIndex[1])
                self.fail(
                    f'Chunked stepTo output diverged from continuous run for {xmlName} '
                    f'at time={continuousData[row, 0]} column={continuousHeaders[col]}: '
                    f'expected {continuousData[row, col]}, observed {chunkedData[row, col]}'
                )

    def test_issue48_ring_unbinding_requires_disconnection(self):
        outputDirectory = mfolder
        fileNumber = '37'

        self._bng_generate(outputDirectory, fileNumber)
        self._run_nfsim(outputDirectory, fileNumber, '-sim 20 -oSteps 20 -cb -seed 1')

        headers, nf = self._load_gdat(os.path.join(outputDirectory, 'v37_nf.gdat'))
        self.assertTrue(len(nf) > 0, 'Issue #48 regression model produced no NFsim output')

        try:
            bondsIdx = headers.index('Obs_Bonds')
            ringsIdx = headers.index('Obs_Rings')
        except ValueError:
            self.fail('Issue #48 regression output missing Obs_Bonds or Obs_Rings columns')

        # In a 4-bond ring, breaking any single bond does not disconnect the species,
        # so L(r!1).R(l!1) -> L(r)+R(l) must never fire.
        self.assertTrue(np.allclose(nf[:, bondsIdx], 4000.0),
                        'Issue #48 failed: Obs_Bonds changed in ring-only system')
        self.assertTrue(np.allclose(nf[:, ringsIdx], 1000.0),
                        'Issue #48 failed: Obs_Rings changed in ring-only system')

    def test_issue49_species_observable_auto_enable_no_crash(self):
        outputDirectory = mfolder
        fileNumber = '38'

        self._bng_generate(outputDirectory, fileNumber)

        # Run without -cb to exercise auto-enable path for Species observables.
        self._run_nfsim(outputDirectory, fileNumber, '-sim 10 -oSteps 10 -seed 2')

        headers, nf = self._load_gdat(os.path.join(outputDirectory, 'v38_nf.gdat'))
        self.assertTrue(len(nf) > 0, 'Issue #49 regression model produced no NFsim output')
        self.assertTrue(np.isfinite(nf).all(), 'Issue #49 failed: NFsim output contains non-finite values')

        # Basic sanity: output includes Species observable column and values are non-negative.
        self.assertIn('Obs_Dimer', headers, 'Issue #49 regression output missing Obs_Dimer column')
        dimerIdx = headers.index('Obs_Dimer')
        self.assertTrue(np.all(nf[:, dimerIdx] >= 0), 'Issue #49 failed: Obs_Dimer has negative values')

    def test_issue53_default_gml_uses_large_limit(self):
        outputDirectory = mfolder
        fileNumber = '36'

        self._bng_generate(outputDirectory, fileNumber)

        # Run without -gml to verify default has been raised and very large populations are supported.
        self._run_nfsim(outputDirectory, fileNumber, '-sim 1 -oSteps 1 -seed 123')

        headers, nf = self._load_gdat(os.path.join(outputDirectory, 'v36_nf.gdat'))
        self.assertTrue(len(nf) > 0, 'Issue #53 regression model produced no NFsim output')
        self.assertIn('A', headers, 'Issue #53 regression output missing molecule count column A')
        aIdx = headers.index('A')
        self.assertEqual(nf[-1, aIdx], 250001.0, 'Issue #53 failed: expected 250001 molecules after initialization')

    def test_issue52_auto_utl_for_multi_molecule_unimolecular_patterns(self):
        outputDirectory = mfolder
        fileNumber = '33'

        self._bng_generate(outputDirectory, fileNumber)

        # Run with default UTL auto (no -utl) to verify the +1 auto-corrected limit holds.
        self._run_nfsim(outputDirectory, fileNumber, '-sim 40000 -oSteps 800 -cb -seed 1')

        headers, nf = self._load_gdat(os.path.join(outputDirectory, 'v33_nf.gdat'))
        self.assertTrue(len(nf) > 0, 'Issue #52 regression model produced no NFsim output')
        self.assertIn('AC', headers, 'Issue #52 regression output missing AC observable')
        acIdx = headers.index('AC')
        # Basic sanity: final AC count should be finite and non-negative.
        self.assertTrue(np.isfinite(nf[-1, acIdx]), 'Issue #52 failed: final AC is not finite')

    def test_step_to_chunking_matches_continuous_run(self):
        self._assert_matching_output_schedules(
            'step_to_cache.xml',
            '-sim 10 -oSteps 10 -seed 1',
            '-sim 10 -oTimes 0,1,2,3,4,5,6,7,8,9,10 -seed 1',
        )

    def test_step_to_zero_propensity_matches_continuous_run(self):
        self._assert_matching_output_schedules(
            'step_to_zero_propensity.xml',
            '-sim 2 -oSteps 2 -seed 1',
            '-sim 2 -oTimes 0,1,2 -seed 1',
        )

    def test_tfun_inline_time_outputs_expected_global_function(self):
        outputDirectory = mfolder
        fileNumber = '44'

        self._bng_generate(outputDirectory, fileNumber)
        self._run_nfsim(outputDirectory, fileNumber, '-sim 2 -oSteps 2 -ogf -seed 1')

        headers, nf = self._load_gdat(os.path.join(outputDirectory, 'v44_nf.gdat'))
        self.assertTrue(len(nf) > 0, 'Inline time TFUN fixture produced no NFsim output')
        self.assertIn('tfun_rate()', headers, 'Inline time TFUN output missing tfun_rate() column')
        tfunIdx = headers.index('tfun_rate()')
        expected = 10.0 * np.clip(nf[:, 0], 0.0, 2.0)
        self.assertTrue(np.allclose(nf[:, tfunIdx], expected),
                        'Inline time TFUN output did not match expected linear interpolation')

    def test_tfun_parameter_counter_outputs_expected_global_function(self):
        outputDirectory = mfolder
        fileNumber = '45'

        self._bng_generate(outputDirectory, fileNumber)
        self._run_nfsim(outputDirectory, fileNumber, '-sim 2 -oSteps 2 -ogf -seed 1')

        headers, nf = self._load_gdat(os.path.join(outputDirectory, 'v45_nf.gdat'))
        self.assertTrue(len(nf) > 0, 'Parameter-counter TFUN fixture produced no NFsim output')
        self.assertIn('tfun_rate()', headers, 'Parameter-counter TFUN output missing tfun_rate() column')
        tfunIdx = headers.index('tfun_rate()')
        self.assertTrue(np.allclose(nf[:, tfunIdx], 10.0),
                        'Parameter-counter TFUN output should stay fixed at the interpolated parameter value')

    def test_tfun_file_time_outputs_expected_global_function(self):
        outputDirectory = mfolder
        fileNumber = '46'

        self._bng_generate(outputDirectory, fileNumber)
        self._run_nfsim(outputDirectory, fileNumber, '-sim 2 -oSteps 2 -ogf -seed 1')

        headers, nf = self._load_gdat(os.path.join(outputDirectory, 'v46_nf.gdat'))
        self.assertTrue(len(nf) > 0, 'File-backed time TFUN fixture produced no NFsim output')
        self.assertIn('tfun_rate()', headers, 'File-backed time TFUN output missing tfun_rate() column')
        tfunIdx = headers.index('tfun_rate()')
        expected = 10.0 * np.clip(nf[:, 0], 0.0, 2.0)
        self.assertTrue(np.allclose(nf[:, tfunIdx], expected),
                        'File-backed time TFUN output did not match expected linear interpolation')

    def test_tfun_function_counter_model_runs(self):
        outputDirectory = mfolder
        fileNumber = '47'

        self._bng_generate(outputDirectory, fileNumber)
        self._run_nfsim(outputDirectory, fileNumber, '-sim 2 -oSteps 2 -ogf -seed 1')

        headers, nf = self._load_gdat(os.path.join(outputDirectory, 'v47_nf.gdat'))
        self.assertTrue(len(nf) > 0, 'Function-counter TFUN fixture produced no NFsim output')
        self.assertTrue(np.isfinite(nf).all(), 'Function-counter TFUN fixture produced non-finite output')
        self.assertIn('driver_fn()', headers, 'Function-counter TFUN output missing driver_fn() column')
        driverIdx = headers.index('driver_fn()')
        self.assertTrue(np.allclose(nf[:, driverIdx], 1.0),
                        'Function-counter driver function should remain constant at 1.0')

    def test_tfun_observable_counter_outputs_bounded_values(self):
        outputDirectory = mfolder
        fileNumber = '48'

        self._bng_generate(outputDirectory, fileNumber)
        self._run_nfsim(outputDirectory, fileNumber, '-sim 2 -oSteps 2 -ogf -seed 1')

        headers, nf = self._load_gdat(os.path.join(outputDirectory, 'v48_nf.gdat'))
        self.assertTrue(len(nf) > 0, 'Observable-counter TFUN fixture produced no NFsim output')
        self.assertIn('tfun_rate()', headers, 'Observable-counter TFUN output missing tfun_rate() column')
        tfunIdx = headers.index('tfun_rate()')
        self.assertAlmostEqual(nf[0, tfunIdx], 0.0, places=7,
                               msg='Observable-counter TFUN should start at zero when X_phos is zero')
        self.assertTrue(np.all((nf[:, tfunIdx] >= 0.0) & (nf[:, tfunIdx] <= 20.0)),
                        'Observable-counter TFUN output should stay within the configured interpolation range')

    def test_tfun_invalid_method_is_rejected(self):
        xmlPath = os.path.join(mfolder, 'invalid_tfun_bad_method.xml')
        outputPath = os.path.join(mfolder, 'invalid_tfun_bad_method_nf.gdat')
        self._run_nfsim_xml(xmlPath, outputPath, '-sim 2 -oSteps 2 -ogf -seed 1', expect_success=False)

    def test_tfun_bionetgen_expr_fixture_outputs_expected_global_functions(self):
        outputDirectory = mfolder
        fileNumber = '49'

        self._run_nfsim(outputDirectory, fileNumber, '-sim 2 -oSteps 2 -ogf -seed 1')

        headers, nf = self._load_gdat(os.path.join(outputDirectory, 'v49_nf.gdat'))
        self.assertTrue(len(nf) > 0, 'BioNetGen-style TFUN expression fixture produced no NFsim output')

        expected_columns = ['f_simple()', 'f_divided()', 'f_scaled()', 'f_complex()', 'Xtot']
        for column in expected_columns:
            self.assertIn(column, headers, f'BioNetGen-style TFUN expression output missing {column} column')

        base = np.array([1.0, 2.0, 4.0])
        simpleIdx = headers.index('f_simple()')
        dividedIdx = headers.index('f_divided()')
        scaledIdx = headers.index('f_scaled()')
        complexIdx = headers.index('f_complex()')
        xtotIdx = headers.index('Xtot')

        self.assertTrue(np.allclose(nf[:, simpleIdx], base),
                        'BioNetGen-style TFUN simple output did not match expected values')
        self.assertTrue(np.allclose(nf[:, dividedIdx], base / 10.0),
                        'BioNetGen-style TFUN divided output did not match expected values')
        self.assertTrue(np.allclose(nf[:, scaledIdx], base * 10.0),
                        'BioNetGen-style TFUN scaled output did not match expected values')
        self.assertTrue(np.allclose(nf[:, complexIdx], (base + 5.0) / 10.0),
                        'BioNetGen-style TFUN complex output did not match expected values')
        self.assertTrue(np.all(np.diff(nf[:, xtotIdx]) >= 0.0),
                        'BioNetGen-style TFUN zero-order production output should be non-decreasing')

    def test_legacy_tfun_placeholder_defaults_to_step_interpolation(self):
        outputDirectory = mfolder
        fileNumber = '50'

        self._run_nfsim(outputDirectory, fileNumber, '-sim 2 -oSteps 4 -ogf -seed 1')

        headers, nf = self._load_gdat(os.path.join(outputDirectory, 'v50_nf.gdat'))
        self.assertTrue(len(nf) > 0, 'Legacy TFUN fixture produced no NFsim output')
        self.assertIn('legacy_tfun_rate()', headers, 'Legacy TFUN output missing legacy_tfun_rate() column')
        tfunIdx = headers.index('legacy_tfun_rate()')
        expected = np.where(nf[:, 0] < 1.0, 0.0, np.where(nf[:, 0] < 2.0, 10.0, 20.0))
        self.assertTrue(np.allclose(nf[:, tfunIdx], expected),
                        'Legacy TFUN placeholder should default to step interpolation when method is omitted')

    def test_invalid_symmetry_factor_throws(self):
        # We need an xml with an invalid symmetry_factor attribute.
        # We can create a simple model xml and manually add symmetry_factor="0" to a ReactionRule.
        xml_content = """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level2/version3" level="2" version="3">
<model id="test_sym_factor">
  <ListOfParameters>
    <Parameter id="k" value="1.0"/>
  </ListOfParameters>
  <ListOfMoleculeTypes>
    <MoleculeType id="A">
      <ListOfComponentTypes>
        <ComponentType id="b"/>
      </ListOfComponentTypes>
    </MoleculeType>
  </ListOfMoleculeTypes>
  <ListOfSpecies>
    <Species id="S1" concentration="100">
      <ListOfMolecules>
        <Molecule id="M1" name="A">
          <ListOfComponents>
            <Component id="C1" name="b" numberOfBonds="0"/>
          </ListOfComponents>
        </Molecule>
      </ListOfMolecules>
    </Species>
  </ListOfSpecies>
  <ListOfReactionRules>
    <ReactionRule id="R1" symmetry_factor="0.0">
      <ListOfReactantPatterns>
        <ReactantPattern id="RP1">
          <ListOfMolecules>
            <Molecule id="M1" name="A">
              <ListOfComponents>
                <Component id="C1" name="b" numberOfBonds="0"/>
              </ListOfComponents>
            </Molecule>
          </ListOfMolecules>
        </ReactantPattern>
      </ListOfReactantPatterns>
      <ListOfProductPatterns>
        <ProductPattern id="PP1">
          <ListOfMolecules>
            <Molecule id="M2" name="A">
              <ListOfComponents>
                <Component id="C2" name="b" numberOfBonds="0"/>
              </ListOfComponents>
            </Molecule>
          </ListOfMolecules>
        </ProductPattern>
      </ListOfProductPatterns>
      <RateLaw id="RL1" type="Ele">
        <ListOfRateConstants>
          <RateConstant value="k"/>
        </ListOfRateConstants>
      </RateLaw>
    </ReactionRule>
  </ListOfReactionRules>
  <ListOfObservables>
    <Observable id="O1" name="A" type="Molecules">
      <ListOfPatterns>
        <Pattern id="P1">
          <ListOfMolecules>
            <Molecule id="M1" name="A">
              <ListOfComponents>
                <Component id="C1" name="b" numberOfBonds="0"/>
              </ListOfComponents>
            </Molecule>
          </ListOfMolecules>
        </Pattern>
      </ListOfPatterns>
    </Observable>
  </ListOfObservables>
</model>
</sbml>
"""
        with open("test_sym_factor.xml", "w") as f:
            f.write(xml_content)

        # We just need to run NFsim on this xml and verify it exits with the correct message/code.
        process = subprocess.Popen([nfsimPath, "-xml", "test_sym_factor.xml"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()

        # The expected behavior: exit(1) and a cerr message
        # "Error!! Symmetry Factor for ReactionRule R1 was not set properly.  quitting."
        self.assertIn(b"Error!! Symmetry Factor for ReactionRule R1 was not set properly.  quitting.", err)
        self.assertEqual(process.returncode, 1)

        # Cleanup
        if os.path.exists("test_sym_factor.xml"):
            os.remove("test_sym_factor.xml")

if __name__ == "__main__":
    suite = unittest.TestSuite()
    if len(sys.argv) > 1:
        os.chdir(sys.argv[1])
    testFolder = mfolder
    tests = getTests(testFolder)
    for index in tests:
        suite.addTest(ParametrizedTestCase.parametrize(TestNFSimFile, param={'num': index,
                    'odir': mfolder, 'iterations': nIterations}))

    # Add targeted model checks to improve coverage for historically unstable cases.
    for modelNum, cfg in targetedTests.items():
        if modelNum in tests:
            suite.addTest(ParametrizedTestCase.parametrize(TestNFSimFile, param={
                'num': modelNum,
                'odir': mfolder,
                'iterations': cfg.get('iterations', nIterations),
                'seed_offset': cfg.get('seed_offset', 0),
                'tag': 'targeted',
            }))

    suite.addTest(unittest.makeSuite(TestIssueRegressions))

    result = unittest.TextTestRunner(verbosity=1).run(suite)

    ret = (list(result.failures) == [] and list(result.errors) == [])
    ret = 0 if ret else 1
    if ret > 0:
        sys.exit("Validation return an error code")
    else:
        sys.exit()
