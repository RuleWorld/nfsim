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

nIterations = 30
nfsimPrePath = ".."
mfolder = "./basicModels"
targetedTests = {
    # Known noisy models get an extra targeted pass with more attempts.
    "16": {"iterations": 60, "seed_offset": 160000, "tol": 0.5},
    "18": {"iterations": 30, "seed_offset": 100000},
    "19": {"iterations": 30, "seed_offset": 200000},
}
if os.name == "nt":
    nfsimPath = os.path.join(nfsimPrePath, "build", "NFsim.exe")
else:
    nfsimPath = os.path.join(nfsimPrePath, "build", "NFsim")


class ParametrizedTestCase(unittest.TestCase):
    """TestCase classes that want to be parametrized should
    inherit from this class.
    """

    def __init__(self, methodName="runTest", param=None):
        super(ParametrizedTestCase, self).__init__(methodName)
        self.param = param

    @staticmethod
    def parametrize(testcase_klass, param=None):
        """Create a suite containing all tests taken from the given
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
            headers = re.sub(r"\s+", " ", line).split(split)

            for line in dataInput:
                nline = re.sub(r"\s+", " ", line.strip()).split(" ")
                try:
                    timeCourse.append([float(x) for x in nline])
                except:
                    print("++++", nline)
        return headers, np.array(timeCourse)
    except IOError:
        print("no file")
        return [], np.array([])


class TestNFSimFile(ParametrizedTestCase):
    def BNGtrajectoryGeneration(self, outputDirectory, fileNumber):
        bngFileName = os.path.join(outputDirectory, "v{0}.bngl".format(fileNumber))
        bionetgen.run(bngFileName, out=outputDirectory, suppress=True)

    def NFsimtrajectoryGeneration(
        self, outputDirectory, fileNumber, runOptions, seed=None
    ):
        runOptions = [x.strip() for x in runOptions.split(" ") if x.strip()]
        if seed is not None and "-seed" not in runOptions:
            runOptions = runOptions + ["-seed", str(seed)]
        with open(os.devnull, "w") as fnull:
            subprocess.check_call(
                [
                    nfsimPath,
                    "-xml",
                    os.path.join(outputDirectory, "v{0}.xml".format(fileNumber)),
                    "-o",
                    os.path.join(outputDirectory, "v{0}_nf.gdat".format(fileNumber)),
                ]
                + runOptions,
                stdout=fnull,
            )

    def _seed_for_iteration(self, index):
        try:
            modelNum = int(self.param["num"])
        except ValueError:
            modelNum = sum([ord(x) for x in str(self.param["num"])])
        seedOffset = int(self.param.get("seed_offset", 0))
        return seedOffset + (modelNum * 1000) + index + 1

    def loadConfigurationFile(self, outputDirectory, fileNumber):
        with open(
            os.path.join(outputDirectory, "r{0}.txt").format(fileNumber), "r"
        ) as f:
            return f.readlines()

    def test_nfsim(self):
        tol = float(self.param.get("tol", 0.35))
        (modelName, runOptions) = self.loadConfigurationFile(
            self.param["odir"], self.param["num"]
        )
        runTag = self.param.get("tag", "default")
        if runTag == "default":
            print(f"Processing model r{self.param['num']}.txt: {modelName.strip()}")
        else:
            print(
                f"Processing model r{self.param['num']}.txt ({runTag}): {modelName.strip()}"
            )
        # here we decide if this is a NFsim only run or not
        if modelName.startswith("NFSIM ONLY"):
            seed = self._seed_for_iteration(0)
            self.BNGtrajectoryGeneration(self.param["odir"], self.param["num"])
            self.NFsimtrajectoryGeneration(
                self.param["odir"], self.param["num"], runOptions, seed=seed
            )
            nfh, nf = loadResults(
                os.path.join(
                    self.param["odir"], "v{0}_nf.gdat".format(self.param["num"])
                ),
                " ",
            )
            # here we just need to make sure we managed to get here without errors
            # assert len(nf) > 0
            self.assertTrue(len(nf) > 0 if type(nf) is list else nf.size > 0)
        else:
            bad = np.array([1])
            lastSeed = None
            for index in range(self.param["iterations"]):
                seed = self._seed_for_iteration(index)
                lastSeed = seed
                print(f"Iteration {index + 1} (seed={seed})")
                self.BNGtrajectoryGeneration(self.param["odir"], self.param["num"])
                self.NFsimtrajectoryGeneration(
                    self.param["odir"], self.param["num"], runOptions, seed=seed
                )
                odeh, ode = loadResults(
                    os.path.join(
                        self.param["odir"], "v{0}_ode.gdat".format(self.param["num"])
                    ),
                    " ",
                )
                ssah, ssa = loadResults(
                    os.path.join(
                        self.param["odir"], "v{0}_ssa.gdat".format(self.param["num"])
                    ),
                    " ",
                )
                nfh, nf = loadResults(
                    os.path.join(
                        self.param["odir"], "v{0}_nf.gdat".format(self.param["num"])
                    ),
                    " ",
                )

                # square root difference per iteration
                ssaDiff = (
                    pow(sum(pow(ode[:, 1:] - ssa[:, 1:], 2)), 0.5)
                    if len(ode) > 0 and len(ssa) > 0
                    else 0
                )
                nfDiff = (
                    pow(sum(pow(ode[:, 1:] - nf[:, 1:], 2)), 0.5)
                    if len(ode) > 0 and len(nf) > 0
                    else 0
                )
                rdiff = nfDiff - ssaDiff - (tol * ssaDiff)
                bad = np.where(rdiff > 0)[0]
                if bad.size > 0:
                    print(
                        f"Observables {bad + 1} did not pass at seed={seed}. Trying again"
                    )
                else:
                    print("Check passed.")
                    break
            self.assertTrue(
                bad.size == 0,
                f"Model r{self.param['num']} failed after {self.param['iterations']} deterministic seeds; "
                f"last seed={lastSeed}, failing observables={bad + 1}",
            )


def getTests(directory):
    """
    Gets a list of bngl files that could be correctly translated in a given 'directory'
    """
    matches = []
    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, "*txt"):
            matches.append("".join(filename.split(".")[0][1:]))
    return sorted(matches)


class TestIssueRegressions(unittest.TestCase):
    def _load_gdat(self, filePath):
        with open(filePath, "r") as f:
            headerLine = re.sub(r"\s+", " ", f.readline().strip())
        headers = [h for h in headerLine.split(" ") if h and h != "#"]
        data = np.loadtxt(filePath, comments="#")
        if data.ndim == 1:
            data = data.reshape(1, -1)
        # Keep only headers that correspond to numeric columns in data.
        if len(headers) > data.shape[1]:
            headers = headers[-data.shape[1] :]
        return headers, data

    def _bng_generate(self, outputDirectory, fileNumber):
        xmlFileName = os.path.join(outputDirectory, "v{0}.xml".format(fileNumber))
        if os.path.exists(xmlFileName):
            # already generated, no need to rerun BNG
            return
        if bionetgen is None:
            self.fail("bionetgen Python package is required to generate XML fixtures")

        bngFileName = os.path.join(outputDirectory, "v{0}.bngl".format(fileNumber))
        bionetgen.run(bngFileName, out=outputDirectory, suppress=True)

    def _run_nfsim_xml(self, xmlPath, outputPath, runOptions, expect_success=True):
        runOptions = [x.strip() for x in runOptions.split(" ") if x.strip()]
        if os.path.exists(outputPath):
            os.remove(outputPath)
        with open(os.devnull, "w") as fnull:
            result = subprocess.run(
                [nfsimPath, "-xml", xmlPath, "-o", outputPath] + runOptions,
                stdout=fnull,
                stderr=fnull,
            )
        if expect_success:
            self.assertEqual(
                result.returncode, 0, f"NFsim failed for XML fixture {xmlPath}"
            )
            self.assertTrue(
                os.path.exists(outputPath),
                f"NFsim did not create expected output file for {xmlPath}",
            )
        else:
            self.assertTrue(
                result.returncode != 0 or not os.path.exists(outputPath),
                f"NFsim unexpectedly succeeded for XML fixture {xmlPath}",
            )
        return result

    def _run_nfsim(self, outputDirectory, fileNumber, runOptions):
        self._run_nfsim_xml(
            os.path.join(outputDirectory, "v{0}.xml".format(fileNumber)),
            os.path.join(outputDirectory, "v{0}_nf.gdat".format(fileNumber)),
            runOptions,
            expect_success=True,
        )

    def _assert_same_seed_connectivity_parity(self, xmlPath, runOptions, label):
        with tempfile.TemporaryDirectory() as tmpdir:
            offPath = os.path.join(tmpdir, "off.gdat")
            onPath = os.path.join(tmpdir, "on.gdat")

            connectOptions = f"{runOptions} -connect".strip()
            self._run_nfsim_xml(xmlPath, offPath, runOptions)
            self._run_nfsim_xml(xmlPath, onPath, connectOptions)

            offHeaders, offData = self._load_gdat(offPath)
            onHeaders, onData = self._load_gdat(onPath)

            self.assertEqual(
                offHeaders,
                onHeaders,
                f"Connectivity regression changed {label} output columns",
            )
            self.assertEqual(
                offData.shape,
                onData.shape,
                f"Connectivity regression changed {label} output shape",
            )
            self.assertTrue(
                np.array_equal(offData, onData),
                f"Connectivity regression changed the same-seed {label} trajectory",
            )

    def test_connectivity_preserves_seeded_tlbr_trajectory(self):
        self._assert_same_seed_connectivity_parity(
            os.path.join(nfsimPrePath, "test", "tlbr", "tlbr.xml"),
            "-sim 1 -oSteps 100 -seed 1",
            "TLBR",
        )

    def test_connectivity_preserves_seeded_local_function_trajectory(self):
        # testSuite/t3 exercises local-function membership updates on a much
        # smaller model than AN_chemotaxis while still reproducing the
        # master-vs-connect divergence fixed by this branch.
        self._assert_same_seed_connectivity_parity(
            os.path.join(nfsimPrePath, "test", "testSuite", "t3.xml"),
            "-sim 1 -oSteps 20 -seed 1",
            "testSuite t3",
        )

    def _assert_matching_output_schedules(
        self, xmlName, continuousOptions, chunkedOptions
    ):
        xmlPath = os.path.join(mfolder, xmlName)
        with tempfile.TemporaryDirectory(prefix="nfsim_step_to_") as tmpdir:
            continuousOutput = os.path.join(tmpdir, "continuous_nf.gdat")
            chunkedOutput = os.path.join(tmpdir, "chunked_nf.gdat")

            self._run_nfsim_xml(xmlPath, continuousOutput, continuousOptions)
            self._run_nfsim_xml(xmlPath, chunkedOutput, chunkedOptions)

            continuousHeaders, continuousData = self._load_gdat(continuousOutput)
            chunkedHeaders, chunkedData = self._load_gdat(chunkedOutput)

            self.assertEqual(
                chunkedHeaders,
                continuousHeaders,
                f"Chunked stepTo output headers differed from continuous run for {xmlName}",
            )
            self.assertEqual(
                chunkedData.shape,
                continuousData.shape,
                f"Chunked stepTo output shape differed from continuous run for {xmlName}",
            )
            if not np.array_equal(chunkedData, continuousData):
                diffIndex = np.argwhere(chunkedData != continuousData)[0]
                row = int(diffIndex[0])
                col = int(diffIndex[1])
                self.fail(
                    f"Chunked stepTo output diverged from continuous run for {xmlName} "
                    f"at time={continuousData[row, 0]} column={continuousHeaders[col]}: "
                    f"expected {continuousData[row, col]}, observed {chunkedData[row, col]}"
                )

    def test_issue48_ring_unbinding_requires_disconnection(self):
        outputDirectory = mfolder
        fileNumber = "37"

        self._bng_generate(outputDirectory, fileNumber)
        self._run_nfsim(outputDirectory, fileNumber, "-sim 20 -oSteps 20 -cb -seed 1")

        headers, nf = self._load_gdat(os.path.join(outputDirectory, "v37_nf.gdat"))
        self.assertTrue(
            len(nf) > 0, "Issue #48 regression model produced no NFsim output"
        )

        try:
            bondsIdx = headers.index("Obs_Bonds")
            ringsIdx = headers.index("Obs_Rings")
        except ValueError:
            self.fail(
                "Issue #48 regression output missing Obs_Bonds or Obs_Rings columns"
            )

        # In a 4-bond ring, breaking any single bond does not disconnect the species,
        # so L(r!1).R(l!1) -> L(r)+R(l) must never fire.
        self.assertTrue(
            np.allclose(nf[:, bondsIdx], 4000.0),
            "Issue #48 failed: Obs_Bonds changed in ring-only system",
        )
        self.assertTrue(
            np.allclose(nf[:, ringsIdx], 1000.0),
            "Issue #48 failed: Obs_Rings changed in ring-only system",
        )

    def test_issue49_species_observable_auto_enable_no_crash(self):
        outputDirectory = mfolder
        fileNumber = "38"

        self._bng_generate(outputDirectory, fileNumber)

        # Run with -no-bscb to exercise auto-enable path for Species observables:
        # complex bookkeeping is otherwise on by default.
        self._run_nfsim(
            outputDirectory, fileNumber, "-sim 10 -oSteps 10 -seed 2 -no-bscb"
        )

        headers, nf = self._load_gdat(os.path.join(outputDirectory, "v38_nf.gdat"))
        self.assertTrue(
            len(nf) > 0, "Issue #49 regression model produced no NFsim output"
        )
        self.assertTrue(
            np.isfinite(nf).all(),
            "Issue #49 failed: NFsim output contains non-finite values",
        )

        # Basic sanity: output includes Species observable column and values are non-negative.
        self.assertIn(
            "Obs_Dimer", headers, "Issue #49 regression output missing Obs_Dimer column"
        )
        dimerIdx = headers.index("Obs_Dimer")
        self.assertTrue(
            np.all(nf[:, dimerIdx] >= 0),
            "Issue #49 failed: Obs_Dimer has negative values",
        )

    def test_internal_bond_reconnection_runs_with_complex_bookkeeping(self):
        # basicModels/v31 breaks a molecule's internal bond and reconnects one end
        # to another molecule, which can leave a picked MappingSet with an unmapped
        # entry. ReactionClass::fire treats that as a null event, but it used to do
        # so only after checkMolecularity had run, and with complex bookkeeping on
        # that check reads the complex of the unmapped molecule. The model ran with
        # no flags and segfaulted with -cb or -bscb.
        outputDirectory = mfolder
        fileNumber = "31"

        self._bng_generate(outputDirectory, fileNumber)
        for flag in ("-cb", "-bscb"):
            self._run_nfsim(
                outputDirectory,
                fileNumber,
                "-sim 1 -oSteps 1 -seed 1 {0}".format(flag),
            )

    def test_issue53_default_gml_uses_large_limit(self):
        outputDirectory = mfolder
        fileNumber = "36"

        self._bng_generate(outputDirectory, fileNumber)

        # Run without -gml to verify default has been raised and very large populations are supported.
        self._run_nfsim(outputDirectory, fileNumber, "-sim 1 -oSteps 1 -seed 123")

        headers, nf = self._load_gdat(os.path.join(outputDirectory, "v36_nf.gdat"))
        self.assertTrue(
            len(nf) > 0, "Issue #53 regression model produced no NFsim output"
        )
        self.assertIn(
            "A", headers, "Issue #53 regression output missing molecule count column A"
        )
        aIdx = headers.index("A")
        self.assertEqual(
            nf[-1, aIdx],
            250001.0,
            "Issue #53 failed: expected 250001 molecules after initialization",
        )

    def test_issue52_auto_utl_for_multi_molecule_unimolecular_patterns(self):
        outputDirectory = mfolder
        fileNumber = "33"

        self._bng_generate(outputDirectory, fileNumber)

        # Run with default UTL auto (no -utl) to verify the +1 auto-corrected limit holds.
        self._run_nfsim(
            outputDirectory, fileNumber, "-sim 40000 -oSteps 800 -cb -seed 1"
        )

        headers, nf = self._load_gdat(os.path.join(outputDirectory, "v33_nf.gdat"))
        self.assertTrue(
            len(nf) > 0, "Issue #52 regression model produced no NFsim output"
        )
        self.assertIn(
            "AC", headers, "Issue #52 regression output missing AC observable"
        )
        acIdx = headers.index("AC")
        # Basic sanity: final AC count should be finite and non-negative.
        self.assertTrue(
            np.isfinite(nf[-1, acIdx]), "Issue #52 failed: final AC is not finite"
        )

    def test_step_to_chunking_matches_continuous_run(self):
        self._assert_matching_output_schedules(
            "step_to_cache.xml",
            "-sim 10 -oSteps 10 -seed 1",
            "-sim 10 -oTimes 0,1,2,3,4,5,6,7,8,9,10 -seed 1",
        )

    def test_step_to_zero_propensity_matches_continuous_run(self):
        self._assert_matching_output_schedules(
            "step_to_zero_propensity.xml",
            "-sim 2 -oSteps 2 -seed 1",
            "-sim 2 -oTimes 0,1,2 -seed 1",
        )

    def test_mm_rate_law_survives_small_km_with_enzyme_in_excess(self):
        # BioNetGen issue 323. The free substrate of the MM rate law used to be
        # computed as 0.5*(b + sqrt(b*b + 4*Km*S)) with b = S - Km - E. That form
        # cancels catastrophically for b < 0 and rounds to exactly zero once
        # 4*Km*S drops below about 1e-16*b*b, leaving a propensity of zero and a
        # reaction that never fires. With Km << E - S the free substrate is
        # Km*S/(E-S) and the propensity reduces to kcat*S, so the substrate decays
        # as a first order death process: P(t) = S0*(1 - exp(-kcat*t)).
        xmlPath = os.path.join(nfsimPrePath, "test", "MM", "mm_small_km.xml")
        S0, kcat, nSeeds = 100.0, 1.0, 20

        total = None
        headers = None
        with tempfile.TemporaryDirectory() as tmpdir:
            for seed in range(1, nSeeds + 1):
                outPath = os.path.join(tmpdir, "mm_{0}.gdat".format(seed))
                self._run_nfsim_xml(xmlPath, outPath, "-sim 4 -oSteps 4 -seed {0}".format(seed))
                headers, data = self._load_gdat(outPath)
                total = data if total is None else total + data
        mean = total / nSeeds

        times = mean[:, headers.index("time")]
        meanP = mean[:, headers.index("Pn")]
        expected = S0 * (1.0 - np.exp(-kcat * times))

        # The reaction fired at all. On the old expression every seed returns a
        # flat zero trajectory, so this alone catches the regression.
        self.assertTrue(
            meanP[-1] > 0.0,
            "MM rate law never fired for small Km with the enzyme in excess",
        )
        # Well clear of the stochastic spread: at t=1 the standard error of the
        # mean over 20 seeds is about 1.7% of the expected value.
        for t, got, want in zip(times[1:], meanP[1:], expected[1:]):
            self.assertAlmostEqual(
                got / want,
                1.0,
                delta=0.15,
                msg="MM product at t={0} averaged {1} over {2} seeds, expected about {3}".format(
                    t, got, nSeeds, want
                ),
            )

    def _mean_final_row(self, xmlPath, runOptions, nSeeds):
        """Mean of the final output row over seeds 1..nSeeds. Returns (headers, row)."""
        total = None
        headers = None
        with tempfile.TemporaryDirectory() as tmpdir:
            for seed in range(1, nSeeds + 1):
                outPath = os.path.join(tmpdir, "run_{0}.gdat".format(seed))
                self._run_nfsim_xml(
                    xmlPath, outPath, "{0} -seed {1}".format(runOptions, seed)
                )
                headers, data = self._load_gdat(outPath)
                total = data if total is None else total + data
        return headers, (total / nSeeds)[-1]

    def test_cb_enforces_same_complex_reactant_molecularity(self):
        # Issue #51. The XML runner routes BioNetGen's complex=>1 / -cb policy
        # into the same distinct-complex molecularity guard as -bscb, and that
        # guard is on by default. A and B below are already in one complex, so
        # their free x/y sites must remain unreacted with no flag and under
        # either complex flag, and may bind only once -no-bscb turns the guard
        # off.
        xmlPath = os.path.join(nfsimPrePath, "test", "Issue51", "issue51.xml")
        with tempfile.TemporaryDirectory() as tmpdir:
            finalCounts = {}
            for label, options in [
                ("none", ""),
                ("cb", "-cb"),
                ("bscb", "-bscb"),
                ("nobscb", "-no-bscb"),
            ]:
                outputPath = os.path.join(tmpdir, "{0}.gdat".format(label))
                self._run_nfsim_xml(
                    xmlPath,
                    outputPath,
                    "-sim 1 -oSteps 1 -seed 1 {0}".format(options),
                )
                headers, data = self._load_gdat(outputPath)
                finalCounts[label] = data[-1, headers.index("bound")]

        self.assertEqual(
            finalCounts["nobscb"],
            1.0,
            "Issue #51 fixture no longer demonstrates same-complex binding",
        )
        self.assertEqual(
            finalCounts["none"],
            0.0,
            "Issue #51 failed: distinct reactant complexes are not enforced by default",
        )
        self.assertEqual(
            finalCounts["cb"],
            0.0,
            "Issue #51 failed: -cb did not enforce distinct reactant complexes",
        )
        self.assertEqual(
            finalCounts["bscb"],
            0.0,
            "Issue #51 failed: -bscb allowed a same-complex reactant match",
        )

    def test_absolute_start_time_offsets_time_functions_and_output(self):
        # Issue #78. The CLI accepted no usable absolute start time: output was
        # still zero-based, and generic BioNetGen time() expressions were not
        # bound to the simulator clock. Explicit output times exercise both the
        # time axis and function value at the requested sample times.
        xmlPath = os.path.join(nfsimPrePath, "test", "Issue78", "issue78.xml")
        with tempfile.TemporaryDirectory() as tmpdir:
            outputPath = os.path.join(tmpdir, "issue78.gdat")
            self._run_nfsim_xml(
                xmlPath,
                outputPath,
                "-i 100 -sim 2 -oTimes 0,1,2 -ogf -seed 1",
            )
            headers, data = self._load_gdat(outputPath)

        timeIdx = headers.index("time")
        stimulusIdx = headers.index("stimulus()")
        xIdx = headers.index("Xtot")
        self.assertTrue(
            np.allclose(data[:, timeIdx], [100.0, 101.0, 102.0]),
            "Issue #78 failed: output time axis did not include the absolute start time",
        )
        self.assertTrue(
            np.allclose(data[:, stimulusIdx], data[:, timeIdx]),
            "Issue #78 failed: generic time() function did not follow the simulator clock",
        )
        self.assertTrue(
            np.all(np.diff(data[:, xIdx]) >= 0.0),
            "Issue #78 fixture produced decreasing zero-order product counts",
        )

    def test_species_observable_function_dependency_updates_after_events(self):
        # Issue #86. Species observable counts were changed through the bulk
        # straightAdd/straightSubtract path during events, so functional rates
        # depending on them stayed at their initial value until output rebuilt
        # the count. The paired production reactions must therefore agree when
        # one uses the Species count and the other uses the molecule count.
        xmlPath = os.path.join(nfsimPrePath, "test", "Issue86", "issue86.xml")
        headers, mean = self._mean_final_row(
            xmlPath, "-sim 300 -oSteps 1 -cb", 5
        )

        self.assertAlmostEqual(
            mean[headers.index("Sobs")],
            mean[headers.index("Mobs")],
            delta=1e-6,
            msg="Issue #86 fixture's Species and molecule observables diverged",
        )
        speciesProduct = mean[headers.index("Ps_n")]
        moleculeProduct = mean[headers.index("Pm_n")]
        self.assertGreater(moleculeProduct, 0.0)
        self.assertLess(
            speciesProduct / moleculeProduct,
            1.5,
            "Issue #86 failed: Species-observable-dependent production stayed stale",
        )
        self.assertLess(
            abs(speciesProduct - moleculeProduct),
            1000.0,
            "Issue #86 failed: Species-observable-dependent rate did not refresh",
        )

    def test_totalrate_is_rejected_for_unsupported_rate_laws(self):
        # Issue #91. These rate-law implementations have no TotalRate
        # semantics. They must fail explicitly instead of silently running the
        # per-instance interpretation.
        cases = [
            (
                os.path.join(nfsimPrePath, "test", "MM", "mm_small_km.xml"),
                '<RateLaw id="RR1_RateLaw" type="MM" totalrate="0">',
                '<RateLaw id="RR1_RateLaw" type="MM" totalrate="1">',
                "TotalRate keyword is not compatible with MM",
            ),
            (
                os.path.join(nfsimPrePath, "test", "testSuite", "t3.xml"),
                '<RateLaw id="RR9_RateLaw" type="Function" name="_rateLaw1" totalrate="0">',
                '<RateLaw id="RR9_RateLaw" type="Function" name="_rateLaw1" totalrate="1">',
                "TotalRate keyword is not compatible with local functions",
            ),
        ]

        for xmlPath, sourceText, replacement, errorText in cases:
            with tempfile.TemporaryDirectory() as tmpdir:
                invalidPath = os.path.join(tmpdir, "invalid.xml")
                outputPath = os.path.join(tmpdir, "invalid.gdat")
                with open(xmlPath, "r") as source:
                    xml = source.read()
                self.assertIn(sourceText, xml, "TotalRate fixture no longer matches source XML")
                with open(invalidPath, "w") as invalid:
                    invalid.write(xml.replace(sourceText, replacement, 1))

                result = subprocess.run(
                    [nfsimPath, "-xml", invalidPath, "-o", outputPath, "-sim", "1"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True,
                )

                self.assertTrue(
                    result.returncode != 0 or not os.path.exists(outputPath),
                    "unsupported TotalRate input unexpectedly produced a trajectory",
                )
                self.assertIn(
                    errorText,
                    result.stdout + result.stderr,
                    "unsupported TotalRate input did not report a clear error",
                )

    def test_symmetry_factor_is_applied_on_every_rate_law(self):
        # ReactionClass's constructor scaled its own baseRate *argument*, which the
        # member had already been copied out of, so the symmetry correction was
        # computed and discarded. Only Ele recovered it, via setBaseRate(). Every
        # other rate law is built with baseRate=1 and never calls setBaseRate, so a
        # rule whose reactant pattern has a non-trivial automorphism ran at
        # 1/symmetry_factor times its intended rate -- 2x for a homodimer.
        #
        # Five dimer pools of X0 copies decay under five rate laws that all encode
        # the same intended per-dimer rate mu, so one expected value covers them
        # all: X0*exp(-mu*t) = 4000*exp(-1) = 1471.5. A dropped factor gives
        # 4000*exp(-2) = 541.3 instead, which is nowhere near it.
        xmlPath = os.path.join(
            nfsimPrePath, "test", "symmetry", "symmetry_factor_rate_laws.xml"
        )
        headers, mean = self._mean_final_row(xmlPath, "-sim 1000 -oSteps 2 -cb", 5)

        expected = 4000.0 * np.exp(-1.0)
        # Per-seed scatter is ~30 counts, so a 5-seed mean has sigma ~13. A +/-170
        # band is ~13 sigma wide and still leaves the dropped-factor value 930
        # counts outside it.
        tolerance = 170.0

        pools = [
            ("Sym_fn", "symmetric, global function (FunctionalRxnClass)"),
            ("Sym_dor", "symmetric, local function (DORRxnClass)"),
            ("Sym_mm", "symmetric, Michaelis-Menten (MMRxnClass)"),
            # These three were already correct; they guard against applying the
            # factor twice, or applying it to an asymmetric rule.
            ("Sym_k", "symmetric, constant rate (BasicRxnClass, control)"),
            ("Asym_fn", "asymmetric, global function (control)"),
            ("Asym_mm", "asymmetric, Michaelis-Menten (control)"),
        ]
        for name, description in pools:
            got = mean[headers.index(name)]
            self.assertAlmostEqual(
                got,
                expected,
                delta=tolerance,
                msg="{0} ({1}) ended at {2:.1f}, expected about {3:.1f}; "
                "{4:.1f} means the symmetry factor was dropped".format(
                    name, description, got, expected, 4000.0 * np.exp(-2.0)
                ),
            )

        # The sharpest form of the claim, needing no expected value at all: a
        # rule's decay must not depend on whether its reactant pattern happens
        # to be symmetric.
        for symName, asymName, rateLaw in [
            ("Sym_fn", "Asym_fn", "global function"),
            ("Sym_mm", "Asym_mm", "Michaelis-Menten"),
        ]:
            sym = mean[headers.index(symName)]
            asym = mean[headers.index(asymName)]
            self.assertAlmostEqual(
                sym,
                asym,
                delta=2 * tolerance,
                msg="{0}: symmetric pool ended at {1:.1f}, asymmetric at {2:.1f}".format(
                    rateLaw, sym, asym
                ),
            )

    def test_totalrate_rules_are_not_scaled_by_the_symmetry_factor(self):
        # FunctionalRxnClass::update_a() scales the propensity by baseRate, which
        # for that class is the reaction center symmetry factor and nothing else.
        # It applied it whether or not the rule uses TotalRate, so a TotalRate
        # rule with a symmetric reaction center ran at a fraction of the rate the
        # model asks for, one half on a homodimer.
        #
        # TotalRate means the rate law gives the whole propensity of the rule. The
        # symmetry factor corrects a match count, and under TotalRate there is no
        # count to correct, so the factor must not be applied.
        #
        # Every BNG-generated TotalRate rule reaches this class and no other.
        # BNG2.pl forces a TotalRate rate law into a Function even when it is a
        # bare constant, which routes it past the Ele/setBaseRate() path, and it
        # rejects TotalRate on Sat/MM/Hill, on Arrhenius, and on local functions.
        #
        # Under TotalRate the propensity is constant while both reactant lists are
        # non-empty, so consumption is linear in time and the expected survivor
        # count is closed form. Both pools are built to land on the same one:
        # Tsym fires kt*T times consuming two A each, Tasym fires 2*kt*T times
        # consuming one B each, so both leave X0 - 2*kt*T = 4000 - 2000 = 2000.
        # A wrongly applied 0.5 halves Tsym's propensity and leaves 3000 instead.
        xmlPath = os.path.join(
            nfsimPrePath, "test", "symmetry", "symmetry_factor_total_rate.xml"
        )
        headers, mean = self._mean_final_row(xmlPath, "-sim 1000 -oSteps 2", 5)

        expected = 2000.0
        halved = 3000.0
        # Firing counts are Poisson, so the symmetric pool's per-seed scatter is
        # 2*sqrt(1000) ~= 63 counts and a 5-seed mean has sigma ~= 28. A +/-200
        # band is ~7 sigma wide and still leaves the halved value 800 counts
        # outside it.
        tolerance = 200.0

        pools = [
            ("Tsym_free", "symmetric reaction center, TotalRate"),
            # symmetry_factor is 1 here, so no placement of the factor can move
            # this pool. It fails only if TotalRate handling broke some other way.
            ("Tasym_free", "asymmetric reaction center, TotalRate (control)"),
        ]
        for name, description in pools:
            got = mean[headers.index(name)]
            self.assertAlmostEqual(
                got,
                expected,
                delta=tolerance,
                msg="{0} ({1}) ended at {2:.1f}, expected about {3:.1f}; "
                "{4:.1f} means the symmetry factor was applied to a rule that "
                "states its own total rate".format(
                    name, description, got, expected, halved
                ),
            )

        # The sharpest form of the claim, needing no expected value at all: under
        # TotalRate a rule's propensity is whatever it says it is, so it must not
        # depend on whether its reactant pattern happens to be symmetric.
        sym = mean[headers.index("Tsym_free")]
        asym = mean[headers.index("Tasym_free")]
        self.assertAlmostEqual(
            sym,
            asym,
            delta=2 * tolerance,
            msg="TotalRate: symmetric pool ended at {0:.1f}, asymmetric at "
            "{1:.1f}".format(sym, asym),
        )

    def test_michaelis_menten_symmetry_factor_scales_the_substrate_count(self):
        # Where the MM law is linear in the substrate match count, scaling that
        # count and scaling the finished propensity coincide, so the fixture above
        # cannot tell them apart. This one sits at X0/Km = 0.4, where they
        # separate. The two rules differ only in whether the substrate dimer's
        # halves are the same molecule type, so they must decay together; there is
        # no closed form here and the pairing is the oracle.
        #
        # Measured over 10 seeds at t=2000:
        #     no factor at all              sym 164.4   asym 761.2
        #     factor on the propensity      sym 993.9   asym 751.8
        #     factor on the substrate count sym 758.0   asym 756.3
        xmlPath = os.path.join(
            nfsimPrePath, "test", "symmetry", "symmetry_factor_mm_saturated.xml"
        )
        headers, mean = self._mean_final_row(xmlPath, "-sim 2000 -oSteps 2 -cb", 10)

        sym = mean[headers.index("Sym_mm")]
        asym = mean[headers.index("Asym_mm")]
        # Per-seed scatter is ~25 counts, so a 10-seed mean has sigma ~8. An
        # 80-count band is ~10 sigma wide and leaves the propensity placement's
        # ~240 gap far outside.
        self.assertAlmostEqual(
            sym,
            asym,
            delta=80.0,
            msg="saturated MM: symmetric ended at {0:.1f}, asymmetric at {1:.1f}. A "
            "gap near +240 means the factor is being applied to the finished "
            "propensity instead of to the substrate count".format(sym, asym),
        )
        # Guard the fixture itself: if a parameter edit drifted this model back
        # into the linear regime the assertion above would keep passing while
        # having stopped discriminating. Linear-regime decay leaves ~1471.
        self.assertLess(
            asym,
            1100.0,
            "the MM control ended at {0:.1f}, too close to the linear-regime value "
            "-- this fixture no longer probes the nonlinear range".format(asym),
        )

    def test_multibond_ring_opening_dissociation_can_fire(self):
        # Product molecularity for a unimolecular unbinding rule used to be tested
        # one deleted bond at a time. For a rule that opens a cyclic complex by
        # deleting several bonds at once that is the wrong question: each ring bond
        # alone leaves the partners connected through the others, so the check
        # refused every one of them and the dissociation never fired.
        #
        # 197 copies of the two-bond ring M(h!1,f!2).M(h!2,f!1), whose only reaction
        # is the reverse homodimerization deleting both ring bonds. BNG's
        # generate_network() integrated as ODEs relaxes to about 63 free monomers.
        # Before the fix the monomer count stayed pinned at 0.
        xmlPath = os.path.join(
            nfsimPrePath, "test", "molecularity", "ring2_homodimer.xml"
        )
        headers, mean = self._mean_final_row(xmlPath, "-sim 200000 -oSteps 2 -bscb", 3)

        # Conservation first: 197 rings is 394 monomer-equivalents of M.
        self.assertAlmostEqual(mean[headers.index("Mtot")], 394.0, delta=1e-6)

        monomers = mean[headers.index("monomers")]
        # The decisive contrast is 0 (trapped) against the ~63 equilibrium, so a
        # generous band still separates the two hypotheses completely.
        self.assertGreater(
            monomers,
            30.0,
            "the two-bond ring did not dissociate (monomers={0:.1f}); the "
            "product-molecularity check is still being applied one bond at a "
            "time".format(monomers),
        )
        self.assertLess(monomers, 110.0)

    def test_singlebond_ring_dissociation_stays_blocked(self):
        # The negative control, and the behavior issues #48 and #61 produced:
        # 100 copies of a size-2 ring whose rule deletes a single L-R bond, which
        # leaves the partners connected through the rest of the cycle. The products
        # do not separate, the network generator drops the reaction, and the ring
        # count must stay at 100 exactly -- no variance across seeds.
        xmlPath = os.path.join(
            nfsimPrePath, "test", "molecularity", "ring_singlebond.xml"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            for seed in (1, 2, 3):
                outPath = os.path.join(tmpdir, "ring_{0}.gdat".format(seed))
                self._run_nfsim_xml(
                    xmlPath, outPath, "-sim 5 -oSteps 1 -bscb -seed {0}".format(seed)
                )
                headers, data = self._load_gdat(outPath)
                rings = data[-1][headers.index("Ring2")]
                self.assertAlmostEqual(
                    rings,
                    100.0,
                    delta=1e-6,
                    msg="a single-bond break inside a ring fired (Ring2={0}, "
                    "seed={1}); product-molecularity enforcement has "
                    "regressed".format(rings, seed),
                )

    def test_ring_dissociation_stays_blocked_when_the_map_has_no_targets(self):
        # NFsim decides whether a deleted bond must separate its two molecules by
        # following them through the rule's Map into the product patterns. XML
        # written by old versions of BioNetGen (test/tlbr/tlbr.xml is one) gives
        # its MapItems no targetID, so there is nothing to follow. NFsim must then
        # fall back on the rule it used before: one reactant pattern and several
        # product patterns means every deleted bond separates its molecules.
        #
        # This strips the targets from ring_singlebond.xml and requires the ring
        # count to stay at 100, as it does with the targets in place.
        sourcePath = os.path.join(
            nfsimPrePath, "test", "molecularity", "ring_singlebond.xml"
        )
        with open(sourcePath, "r") as source:
            xml = source.read()
        stripped, nTargets = re.subn(r' targetID="[^"]*"', "", xml)
        self.assertGreater(nTargets, 0, "ring_singlebond.xml has no map targets")
        with tempfile.TemporaryDirectory() as tmpdir:
            xmlPath = os.path.join(tmpdir, "ring_singlebond_no_targets.xml")
            with open(xmlPath, "w") as target:
                target.write(stripped)
            for seed in (1, 2, 3):
                outPath = os.path.join(tmpdir, "ring_{0}.gdat".format(seed))
                self._run_nfsim_xml(
                    xmlPath, outPath, "-sim 5 -oSteps 1 -seed {0}".format(seed)
                )
                headers, data = self._load_gdat(outPath)
                rings = data[-1][headers.index("Ring2")]
                self.assertAlmostEqual(
                    rings,
                    100.0,
                    delta=1e-6,
                    msg="a ring bond was broken when the rule's map had no targets "
                    "(Ring2={0}, seed={1})".format(rings, seed),
                )

    def test_species_observable_counts_without_a_bookkeeping_flag(self):
        # System.useComplex used to be derived solely from blockSameComplexBinding.
        # A Species-typed observable is tallied by iterating complexes, so run
        # without complex bookkeeping it was counted with tracking disabled and
        # reported a number with no physical ceiling.
        #
        # ring2_homodimer.xml declares the Species observable `dimers` over 197
        # two-bond rings, so the count can never exceed 197. Run with no
        # bookkeeping flag, it used to report ~3612. It must now be physical, and
        # agree with the bookkeeping-enabled run -- this model has no same-complex
        # binding, so the two modes describe the same trajectory.
        #
        # Bookkeeping is on by default now, so the run in which nothing but the
        # observable asks for it is the one given -no-bscb.
        xmlPath = os.path.join(
            nfsimPrePath, "test", "molecularity", "ring2_homodimer.xml"
        )
        headers, checksOff = self._mean_final_row(
            xmlPath, "-sim 200000 -oSteps 2 -no-bscb", 3
        )
        _, withBscb = self._mean_final_row(xmlPath, "-sim 200000 -oSteps 2 -bscb", 3)

        seeded = 197.0
        dimers = checksOff[headers.index("dimers")]
        self.assertLessEqual(
            dimers,
            seeded,
            "the Species observable reported {0:.1f} dimers with bookkeeping not "
            "requested, but only {1:.0f} rings were ever seeded -- complex tracking "
            "is not on for a model that needs it".format(dimers, seeded),
        )
        self.assertAlmostEqual(
            dimers,
            withBscb[headers.index("dimers")],
            delta=20.0,
            msg="the Species observable disagrees between the -no-bscb run and the "
            "-bscb run, which describe the same trajectory for this model",
        )
        # The molecules themselves were always counted correctly; it is only the
        # complex-level tally that was wrong. Guards against "fixing" the count by
        # perturbing the underlying simulation.
        self.assertAlmostEqual(checksOff[headers.index("Mtot")], 394.0, delta=1e-6)

    def test_ring_bond_dissociation_is_blocked_by_default(self):
        # lanl/bngsim#632. "A.B -> A + B" may fire only when the bonds it deletes
        # separate A from B, and "A + B -> A.B" only when A and B are in different
        # complexes. NFsim applied both checks only when -cb or -bscb was given, so
        # a model run with neither flag broke ring bonds and reported plausible,
        # wrong counts with no warning. Both are now on by default, as they are in
        # BioNetGen's network generator, and -no-bscb turns them off together.
        #
        # ring_closure.xml holds 500 bivalent L and 1000 bivalent R. An L-R pair can
        # close into a two-bond ring and no rule opens one, so the ring is absorbing
        # and every L must end in it: Ring is 500 exactly, at every seed. Ten seeds
        # all reach 500 by t=100. With the checks off, rings open as fast as they
        # close and Ring settles near 50 (48 and 52 at seeds 1 and 2).
        xmlPath = os.path.join(nfsimPrePath, "test", "molecularity", "ring_closure.xml")
        with tempfile.TemporaryDirectory() as tmpdir:
            for label, options in [("default", ""), ("cb", "-cb"), ("bscb", "-bscb")]:
                for seed in (1, 2, 3):
                    outPath = os.path.join(tmpdir, "{0}_{1}.gdat".format(label, seed))
                    self._run_nfsim_xml(
                        xmlPath,
                        outPath,
                        "-sim 1000 -oSteps 2 -seed {0} {1}".format(seed, options),
                    )
                    headers, data = self._load_gdat(outPath)
                    rings = data[-1][headers.index("Ring")]
                    self.assertAlmostEqual(
                        rings,
                        500.0,
                        delta=1e-6,
                        msg="Ring ended at {0} instead of 500 ({1}, seed={2}); a "
                        "'+' rule is breaking ring bonds or closing rings inside "
                        "one complex".format(rings, label, seed),
                    )
                    self.assertAlmostEqual(
                        data[-1][headers.index("Lfree")], 0.0, delta=1e-6
                    )

            # The opt-out has to stay reachable, and has to turn off the same
            # thing: this is the count the model gave with no flag before.
            outPath = os.path.join(tmpdir, "checks_off.gdat")
            self._run_nfsim_xml(
                xmlPath, outPath, "-sim 1000 -oSteps 2 -seed 1 -no-bscb"
            )
            headers, data = self._load_gdat(outPath)
            self.assertLess(
                data[-1][headers.index("Ring")],
                250.0,
                "-no-bscb left the molecularity checks on",
            )

    def test_no_bscb_cannot_be_combined_with_cb_or_bscb(self):
        # -no-bscb and -cb/-bscb ask for opposite things. BioNetGen passes -cb
        # unless told complex=>0, so a silent winner in either direction would
        # hand someone the mode they did not ask for.
        xmlPath = os.path.join(nfsimPrePath, "test", "molecularity", "ring_closure.xml")
        with tempfile.TemporaryDirectory() as tmpdir:
            for flag in ("-cb", "-bscb"):
                outPath = os.path.join(tmpdir, "conflict.gdat")
                result = self._run_nfsim_xml(
                    xmlPath,
                    outPath,
                    "-sim 1 -oSteps 1 -seed 1 -no-bscb {0}".format(flag),
                    expect_success=False,
                )
                self.assertEqual(
                    result.returncode,
                    1,
                    "-no-bscb {0} was accepted".format(flag),
                )

    def test_complex_scoped_function_sees_molecules_bound_during_the_run(self):
        # NFsim's reader used to switch complex bookkeeping on at the point where it
        # met a pattern holding two molecules with no bond between them: a disjoint
        # observable or reactant pattern, or the product of any rule that joins its
        # products with ".". By then the seed molecules existed. Each got a complex
        # of its own and kept its cached "not tracked" flag, so those complexes were
        # never merged when the molecules bound, while the System went on evaluating
        # complex-scoped local functions one complex at a time. Adding such a pattern
        # to a model changed its kinetics.
        #
        # Both fixtures bind one B to each A and flip A at k times the number of B in
        # A's complex, so no A flips before it binds and all 100 flip after. With the
        # late switch all 100 bound and none flipped. It takes -no-bscb to show,
        # because bookkeeping is otherwise on from the start.
        for name in ("complex_scope_dot_product", "complex_scope_disjoint_observable"):
            xmlPath = os.path.join(
                nfsimPrePath, "test", "localFunction", "{0}.xml".format(name)
            )
            with tempfile.TemporaryDirectory() as tmpdir:
                for label, options in [("checks_off", "-no-bscb"), ("default", "")]:
                    for seed in (1, 2, 3):
                        outPath = os.path.join(
                            tmpdir, "{0}_{1}.gdat".format(label, seed)
                        )
                        self._run_nfsim_xml(
                            xmlPath,
                            outPath,
                            "-sim 40 -oSteps 2 -seed {0} {1}".format(seed, options),
                        )
                        headers, data = self._load_gdat(outPath)
                        final = dict(zip(headers, data[-1]))
                        self.assertAlmostEqual(final["Abound"], 100.0, delta=1e-6)
                        self.assertAlmostEqual(
                            final["Aon"],
                            100.0,
                            delta=1e-6,
                            msg="{0} ({1}, seed={2}): {3:.0f} of 100 bound A flipped; "
                            "the local function is not seeing the B that A "
                            "bound".format(name, label, seed, final["Aon"]),
                        )

    def test_assistance_neighborhood_model_adapts_with_checks_on(self):
        # an2.bngl holds an enzyme to a receptor cluster by two bonds at once, its
        # tether and its active site. The rules that break one of them used to write
        # "+" between the products although the other bond still joined them. With
        # the molecularity checks on, such a rule can never fire, so the model froze:
        # all 7030 receptor dimers stayed at their seeded methylation level, m~3.
        #
        # Each rule is now split in two, "+" when the other bond is free and "."
        # when it is bound. Ten seeds leave between 6077 and 6242 dimers at m~3 by
        # t=10. Written that way the model no longer depends on the checks, so
        # -no-bscb must give the same trajectory.
        xmlPath = os.path.join(nfsimPrePath, "test", "AN_chemotaxis", "an2.xml")
        levels = ["R{0}".format(m) for m in range(9)]
        with tempfile.TemporaryDirectory() as tmpdir:
            for seed in (1, 2, 3):
                outPath = os.path.join(tmpdir, "an2_{0}.gdat".format(seed))
                self._run_nfsim_xml(
                    xmlPath, outPath, "-sim 10 -oSteps 2 -seed {0}".format(seed)
                )
                headers, data = self._load_gdat(outPath)
                final = dict(zip(headers, data[-1]))
                self.assertAlmostEqual(
                    sum(final[name] for name in levels), 7030.0, delta=1e-6
                )
                self.assertLess(
                    final["R3"],
                    6600.0,
                    "an2 did not adapt (R3={0}, seed={1}); a bond-breaking rule is "
                    "written with '+' while its products stay connected".format(
                        final["R3"], seed
                    ),
                )
                self.assertGreater(final["R3"], 5800.0)

                offPath = os.path.join(tmpdir, "an2_off_{0}.gdat".format(seed))
                self._run_nfsim_xml(
                    xmlPath,
                    offPath,
                    "-sim 10 -oSteps 2 -seed {0} -no-bscb".format(seed),
                )
                _, offData = self._load_gdat(offPath)
                self.assertTrue(
                    np.array_equal(data, offData),
                    "an2 gave a different trajectory with -no-bscb (seed={0})".format(
                        seed
                    ),
                )

    def test_two_reactant_rule_does_not_break_a_ring_bond(self):
        # "A.B + C -> A + B + C" writes "+" between A and B, so it may fire only
        # if deleting the bond leaves them in different complexes. NFsim applied
        # that test only to rules with one reactant pattern, so a rule with two
        # broke ring bonds freely.
        #
        # ring_two_reactants.xml holds 100 A.B pairs joined twice (Ring) and 100
        # joined once (Single). The rule must never fire on a Ring pair and must
        # clear the Single pairs by t=40. BioNetGen's network generator produces
        # exactly that one reaction.
        xmlPath = os.path.join(
            nfsimPrePath, "test", "molecularity", "ring_two_reactants.xml"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            for seed in (1, 2, 3):
                outPath = os.path.join(tmpdir, "two_{0}.gdat".format(seed))
                self._run_nfsim_xml(
                    xmlPath, outPath, "-sim 40 -oSteps 4 -seed {0}".format(seed)
                )
                headers, data = self._load_gdat(outPath)
                rings = data[:, headers.index("Ring")]
                self.assertTrue(
                    np.all(rings == 100.0),
                    "a rule with two reactant patterns broke a ring bond (Ring "
                    "fell to {0}, seed={1})".format(rings.min(), seed),
                )
                self.assertAlmostEqual(
                    data[-1][headers.index("Single")], 0.0, delta=1e-6
                )

    def test_pure_context_reactant_is_counted_once_per_complex(self):
        # Issue #87. BioNetGen gives a reactant pattern the rule does not transform
        # one reaction instance per matching complex, however many molecules inside
        # that complex match it, because every embedding yields the identical
        # reaction. NFsim enumerated matches per molecule, so a homodimeric
        # catalyst drove its rule twice as fast as a heterodimeric one and a
        # homotrimer ring three times as fast.
        #
        # Every pool below carries the same intended per-substrate rate
        # mu = kcat*E0 = 2.5e-4, so the oracle is BNG's own generated network
        # integrated as ODEs: X0*exp(-mu*t) = 4000*exp(-0.5) = 2426.1. The
        # Michaelis-Menten pair has no closed form and lands at 2344.8.
        xmlPath = os.path.join(nfsimPrePath, "test", "context", "context_symmetry.xml")
        headers, mean = self._mean_final_row(xmlPath, "-sim 2000 -oSteps 2 -cb", 10)

        expected = 4000.0 * np.exp(-0.5)
        expectedMM = 2344.8
        # Per-seed scatter is ~30 counts, so a 10-seed mean has sigma ~10. A +/-100
        # band is ~10 sigma wide and leaves every over-counted value (~1471 at 2x,
        # ~893 at 3x) more than 900 counts outside it.
        tolerance = 100.0

        pools = [
            ("Dim_sym", "homodimer catalyst, constant rate", expected),
            ("Fn_sym", "homodimer catalyst, global function", expected),
            ("Dor_sym", "homodimer catalyst, local function", expected),
            ("Mm_sym", "symmetric enzyme, Michaelis-Menten", expectedMM),
            ("Ring_sym", "homotrimer ring catalyst (3x, not 2x)", expected),
            ("Sub_subunit", "single-subunit pattern against a homodimer", expected),
            # The load-bearing case: a single-molecule pattern against a complex
            # holding two *distinguishable* copies. No automorphism anywhere, and
            # still over-counted -- which rules out embeddings/|Aut(pattern)|.
            ("Sub_scaffold", "two distinguishable copies, no symmetry at all", expected),
            # Sharper still: single-molecule pattern against a homotrimer ring, so
            # embeddings/|Aut| would say 3, a hardcoded factor of two 1.5, and
            # once-per-complex 1. BNG says 1.
            ("Sub_trimer", "single-molecule pattern against a homotrimer ring", expected),
        ]
        for name, description, want in pools:
            got = mean[headers.index(name)]
            self.assertAlmostEqual(
                got,
                want,
                delta=tolerance,
                msg="{0} ({1}) ended at {2:.1f}, expected about {3:.1f}; a value "
                "near 1471 means two matching molecules in one complex were "
                "counted as two reaction instances, near 893 means three".format(
                    name, description, got, want
                ),
            )

        # Needs no oracle: two catalysts present at the same complex count and
        # carrying the same rate constant must give the same rate, whatever the
        # rate law.
        for multi, single, rateLaw in [
            ("Dim_sym", "Dim_asym", "constant rate"),
            ("Ring_sym", "Ring_asym", "constant rate, trimer ring"),
            ("Fn_sym", "Fn_asym", "global function"),
            ("Dor_sym", "Dor_asym", "local function"),
            ("Mm_sym", "Mm_asym", "Michaelis-Menten"),
        ]:
            many = mean[headers.index(multi)]
            one = mean[headers.index(single)]
            self.assertAlmostEqual(
                many,
                one,
                delta=2 * tolerance,
                msg="{0}: multi-subunit catalyst ended at {1:.1f}, single-subunit "
                "control at {2:.1f}".format(rateLaw, many, one),
            )

    def test_a_transformed_symmetric_pattern_keeps_both_of_its_sites(self):
        # The other direction, and the reason "which reactants are pure context"
        # has to be decided before finalize() appends its placeholder. Bind_sym's
        # catalyst dimer IS transformed -- one half of it binds -- so its two
        # halves are two genuinely distinct reactive sites and two distinct
        # reactions. BNG agrees explicitly: the generated network gives this rule
        # 2*kb where the heterodimer control gets kb.
        #
        # NFsim marks the second partner of a binding with an EMPTY transform, the
        # same type finalize() uses for its placeholder, so deciding pure context
        # from the transformation types afterwards misclassifies exactly this rule
        # and erases its factor of two -- landing it on top of its own control.
        xmlPath = os.path.join(nfsimPrePath, "test", "context", "context_symmetry.xml")
        headers, mean = self._mean_final_row(xmlPath, "-sim 2000 -oSteps 2 -cb", 10)

        sym = mean[headers.index("Bind_sym")]
        asym = mean[headers.index("Bind_asym")]
        self.assertAlmostEqual(
            sym,
            218.0,
            delta=25.0,
            msg="Bind_sym ended at {0:.1f}, expected about 218. A value near its "
            "control ({1:.1f}) means per-complex counting reached a reactant the "
            "rule transforms and divided out a factor BNG put there on "
            "purpose".format(sym, asym),
        )
        self.assertAlmostEqual(
            asym,
            320.0,
            delta=30.0,
            msg="Bind_asym ended at {0:.1f}, expected about 320. This rule has one "
            "reactive site and one matching molecule per complex, so nothing here "
            "should move it".format(asym),
        )
        # Guard the pair: the two assertions above only discriminate while the
        # two-site arm actually runs faster than its control.
        self.assertGreater(
            asym - sym,
            50.0,
            "Bind_sym and Bind_asym have converged -- this pair no longer tells a "
            "transformed pattern's real multiplicity from a context over-count",
        )

    def test_tfun_inline_time_outputs_expected_global_function(self):
        outputDirectory = mfolder
        fileNumber = "44"

        self._bng_generate(outputDirectory, fileNumber)
        self._run_nfsim(outputDirectory, fileNumber, "-sim 2 -oSteps 2 -ogf -seed 1")

        headers, nf = self._load_gdat(os.path.join(outputDirectory, "v44_nf.gdat"))
        self.assertTrue(
            len(nf) > 0, "Inline time TFUN fixture produced no NFsim output"
        )
        self.assertIn(
            "tfun_rate()", headers, "Inline time TFUN output missing tfun_rate() column"
        )
        tfunIdx = headers.index("tfun_rate()")
        expected = 10.0 * np.clip(nf[:, 0], 0.0, 2.0)
        self.assertTrue(
            np.allclose(nf[:, tfunIdx], expected),
            "Inline time TFUN output did not match expected linear interpolation",
        )

    def test_tfun_parameter_counter_outputs_expected_global_function(self):
        outputDirectory = mfolder
        fileNumber = "45"

        self._bng_generate(outputDirectory, fileNumber)
        self._run_nfsim(outputDirectory, fileNumber, "-sim 2 -oSteps 2 -ogf -seed 1")

        headers, nf = self._load_gdat(os.path.join(outputDirectory, "v45_nf.gdat"))
        self.assertTrue(
            len(nf) > 0, "Parameter-counter TFUN fixture produced no NFsim output"
        )
        self.assertIn(
            "tfun_rate()",
            headers,
            "Parameter-counter TFUN output missing tfun_rate() column",
        )
        tfunIdx = headers.index("tfun_rate()")
        self.assertTrue(
            np.allclose(nf[:, tfunIdx], 10.0),
            "Parameter-counter TFUN output should stay fixed at the interpolated parameter value",
        )

    def test_tfun_file_time_outputs_expected_global_function(self):
        outputDirectory = mfolder
        fileNumber = "46"

        self._bng_generate(outputDirectory, fileNumber)
        self._run_nfsim(outputDirectory, fileNumber, "-sim 2 -oSteps 2 -ogf -seed 1")

        headers, nf = self._load_gdat(os.path.join(outputDirectory, "v46_nf.gdat"))
        self.assertTrue(
            len(nf) > 0, "File-backed time TFUN fixture produced no NFsim output"
        )
        self.assertIn(
            "tfun_rate()",
            headers,
            "File-backed time TFUN output missing tfun_rate() column",
        )
        tfunIdx = headers.index("tfun_rate()")
        expected = 10.0 * np.clip(nf[:, 0], 0.0, 2.0)
        self.assertTrue(
            np.allclose(nf[:, tfunIdx], expected),
            "File-backed time TFUN output did not match expected linear interpolation",
        )

    def test_tfun_function_counter_model_runs(self):
        outputDirectory = mfolder
        fileNumber = "47"

        self._bng_generate(outputDirectory, fileNumber)
        self._run_nfsim(outputDirectory, fileNumber, "-sim 2 -oSteps 2 -ogf -seed 1")

        headers, nf = self._load_gdat(os.path.join(outputDirectory, "v47_nf.gdat"))
        self.assertTrue(
            len(nf) > 0, "Function-counter TFUN fixture produced no NFsim output"
        )
        self.assertTrue(
            np.isfinite(nf).all(),
            "Function-counter TFUN fixture produced non-finite output",
        )
        self.assertIn(
            "driver_fn()",
            headers,
            "Function-counter TFUN output missing driver_fn() column",
        )
        driverIdx = headers.index("driver_fn()")
        self.assertTrue(
            np.allclose(nf[:, driverIdx], 1.0),
            "Function-counter driver function should remain constant at 1.0",
        )

    def test_tfun_observable_counter_outputs_bounded_values(self):
        outputDirectory = mfolder
        fileNumber = "48"

        self._bng_generate(outputDirectory, fileNumber)
        self._run_nfsim(outputDirectory, fileNumber, "-sim 2 -oSteps 2 -ogf -seed 1")

        headers, nf = self._load_gdat(os.path.join(outputDirectory, "v48_nf.gdat"))
        self.assertTrue(
            len(nf) > 0, "Observable-counter TFUN fixture produced no NFsim output"
        )
        self.assertIn(
            "tfun_rate()",
            headers,
            "Observable-counter TFUN output missing tfun_rate() column",
        )
        tfunIdx = headers.index("tfun_rate()")
        self.assertAlmostEqual(
            nf[0, tfunIdx],
            0.0,
            places=7,
            msg="Observable-counter TFUN should start at zero when X_phos is zero",
        )
        self.assertTrue(
            np.all((nf[:, tfunIdx] >= 0.0) & (nf[:, tfunIdx] <= 20.0)),
            "Observable-counter TFUN output should stay within the configured interpolation range",
        )

    def test_tfun_invalid_method_is_rejected(self):
        xmlPath = os.path.join(mfolder, "invalid_tfun_bad_method.xml")
        outputPath = os.path.join(mfolder, "invalid_tfun_bad_method_nf.gdat")
        self._run_nfsim_xml(
            xmlPath, outputPath, "-sim 2 -oSteps 2 -ogf -seed 1", expect_success=False
        )

    def test_tfun_bionetgen_expr_fixture_outputs_expected_global_functions(self):
        outputDirectory = mfolder
        fileNumber = "49"

        self._run_nfsim(outputDirectory, fileNumber, "-sim 2 -oSteps 2 -ogf -seed 1")

        headers, nf = self._load_gdat(os.path.join(outputDirectory, "v49_nf.gdat"))
        self.assertTrue(
            len(nf) > 0,
            "BioNetGen-style TFUN expression fixture produced no NFsim output",
        )

        expected_columns = [
            "f_simple()",
            "f_divided()",
            "f_scaled()",
            "f_complex()",
            "Xtot",
        ]
        for column in expected_columns:
            self.assertIn(
                column,
                headers,
                f"BioNetGen-style TFUN expression output missing {column} column",
            )

        base = np.array([1.0, 2.0, 4.0])
        simpleIdx = headers.index("f_simple()")
        dividedIdx = headers.index("f_divided()")
        scaledIdx = headers.index("f_scaled()")
        complexIdx = headers.index("f_complex()")
        xtotIdx = headers.index("Xtot")

        self.assertTrue(
            np.allclose(nf[:, simpleIdx], base),
            "BioNetGen-style TFUN simple output did not match expected values",
        )
        self.assertTrue(
            np.allclose(nf[:, dividedIdx], base / 10.0),
            "BioNetGen-style TFUN divided output did not match expected values",
        )
        self.assertTrue(
            np.allclose(nf[:, scaledIdx], base * 10.0),
            "BioNetGen-style TFUN scaled output did not match expected values",
        )
        self.assertTrue(
            np.allclose(nf[:, complexIdx], (base + 5.0) / 10.0),
            "BioNetGen-style TFUN complex output did not match expected values",
        )
        self.assertTrue(
            np.all(np.diff(nf[:, xtotIdx]) >= 0.0),
            "BioNetGen-style TFUN zero-order production output should be non-decreasing",
        )

    def test_legacy_tfun_placeholder_defaults_to_step_interpolation(self):
        outputDirectory = mfolder
        fileNumber = "50"

        self._run_nfsim(outputDirectory, fileNumber, "-sim 2 -oSteps 4 -ogf -seed 1")

        headers, nf = self._load_gdat(os.path.join(outputDirectory, "v50_nf.gdat"))
        self.assertTrue(len(nf) > 0, "Legacy TFUN fixture produced no NFsim output")
        self.assertIn(
            "legacy_tfun_rate()",
            headers,
            "Legacy TFUN output missing legacy_tfun_rate() column",
        )
        tfunIdx = headers.index("legacy_tfun_rate()")
        expected = np.where(nf[:, 0] < 1.0, 0.0, np.where(nf[:, 0] < 2.0, 10.0, 20.0))
        self.assertTrue(
            np.allclose(nf[:, tfunIdx], expected),
            "Legacy TFUN placeholder should default to step interpolation when method is omitted",
        )

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
        process = subprocess.Popen(
            [nfsimPath, "-xml", "test_sym_factor.xml"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        out, err = process.communicate()

        # The expected behavior: exit(1) and a cerr message
        # "Error!! Symmetry Factor for ReactionRule R1 was not set properly.  quitting."
        self.assertIn(
            b"Error!! Symmetry Factor for ReactionRule R1 was not set properly.  quitting.",
            err,
        )
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
        suite.addTest(
            ParametrizedTestCase.parametrize(
                TestNFSimFile,
                param={"num": index, "odir": mfolder, "iterations": nIterations},
            )
        )

    # Add targeted model checks to improve coverage for historically unstable cases.
    for modelNum, cfg in targetedTests.items():
        if modelNum in tests:
            suite.addTest(
                ParametrizedTestCase.parametrize(
                    TestNFSimFile,
                    param={
                        "num": modelNum,
                        "odir": mfolder,
                        "iterations": cfg.get("iterations", nIterations),
                        "seed_offset": cfg.get("seed_offset", 0),
                        "tol": cfg.get("tol", 0.35),
                        "tag": "targeted",
                    },
                )
            )

    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestIssueRegressions))

    result = unittest.TextTestRunner(verbosity=1).run(suite)

    ret = list(result.failures) == [] and list(result.errors) == []
    ret = 0 if ret else 1
    if ret > 0:
        sys.exit("Validation return an error code")
    else:
        sys.exit()
