import unittest
import os
import numpy as np
import subprocess
import re
import fnmatch
import sys
import bionetgen

nIterations=15
nfsimPrePath='..'
mfolder='./basicModels'
bngPath = os.path.join(bionetgen.defaults.bng_path, "BNG2.pl")
targetedTests = {
    # Known noisy models get an extra targeted pass with more attempts.
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
            headers = re.sub('\s+', ' ', line).split(split)

            for line in dataInput:
                nline = re.sub('\s+', ' ', line.strip()).split(' ')
                try:
                    timeCourse.append([float(x) for x in nline])
                except:
                    print('++++', nline)
        return headers, np.array(timeCourse)
    except IOError:
        print('no file')
        return [], []


class TestNFSimFile(ParametrizedTestCase):
    # XXX:ideally this should be done through the console but I'm doing the quick and dirty version right now
    def BNGtrajectoryGeneration(self, outputDirectory, fileNumber):
        bngFileName = os.path.join(outputDirectory, 'v{0}.bngl'.format(fileNumber))
        with open(os.devnull, "w") as fnull:
            subprocess.check_call(['perl', bngPath, '-outdir', outputDirectory, '-log', bngFileName], stdout=fnull)

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
            headerLine = re.sub('\s+', ' ', f.readline().strip())
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

        bngFileName = os.path.join(outputDirectory, 'v{0}.bngl'.format(fileNumber))
        with open(os.devnull, "w") as fnull:
            subprocess.check_call(['perl', bngPath, '-outdir', outputDirectory, '-log', bngFileName], stdout=fnull)

    def _run_nfsim(self, outputDirectory, fileNumber, runOptions):
        runOptions = [x.strip() for x in runOptions.split(' ') if x.strip()]
        with open(os.devnull, "w") as fnull:
            subprocess.check_call([
                nfsimPath,
                '-xml', os.path.join(outputDirectory, 'v{0}.xml'.format(fileNumber)),
                '-o', os.path.join(outputDirectory, 'v{0}_nf.gdat'.format(fileNumber))
            ] + runOptions, stdout=fnull)

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
        self.assertGreaterEqual(nf[-1, acIdx], 0.0, 'Issue #52 failed: final AC is negative')

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
