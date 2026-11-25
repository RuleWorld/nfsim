# NFsim - the network free stochastic simulator

[![NFsim build status](https://github.com/RuleWorld/nfsim/workflows/main-validation/badge.svg)](https://github.com/RuleWorld/nfsim/actions)
<a href="https://scan.coverity.com/projects/nfsim">
  <img alt="Coverity Scan Build Status"
       src="https://scan.coverity.com/projects/15734/badge.svg"/>
</a>

## Overview

NFsim is a free, open-source, biochemical reaction simulator designed to handle
systems that have a large or even infinite number of possible molecular
interactions or states. NFsim also has advanced and flexible options for
simulating coarse-grained representations of complex nonlinear reaction
mechanisms.

NFsim is ideal for modeling polymerization, aggregation, and cooperative
reactions that cannot be handled with traditional stochastic or ODE simulators.
Models are specified in the BioNetGen Language, providing a powerful model
building environment.

NFsim is distributed with BioNetGen. To download NFsim on its own, get a
prebuilt binary for Linux, macOS or Windows from the
[releases page](https://github.com/RuleWorld/nfsim/releases). To work on the
code or contribute, fork the repository and open a pull request against
master.

## Citing NFsim

If you use NFsim in your research, please cite both:

1. **This software**, with the version you used: Kutuva A, Hlavacek WS,
   Faeder JR. NFsim, version 2.0.0. https://github.com/RuleWorld/nfsim.
   The "Cite this repository" link on the GitHub page gives this citation in
   APA and BibTeX form (from `CITATION.cff`).
2. **The NFsim paper**, for the network-free simulation method: Sneddon MW,
   Faeder JR, Emonet T. Efficient modeling, simulation and coarse-graining of
   biological complexity with NFsim. *Nature Methods* 8(2):177-183 (2011).
   https://doi.org/10.1038/nmeth.1546

NFsim was created by Michael Sneddon, Justin Hogg, Jose-Juan Tapia, James
Faeder and Thierry Emonet at Yale University and the University of Pittsburgh.
[ACKNOWLEDGMENTS.md](ACKNOWLEDGMENTS.md) credits them and everyone who has
contributed since.

## Repository Contents

NFsim is released under the MIT License. See LICENSE.txt for more details about
redistribution restrictions.

For help with running NFsim, see the user manual, NFsim\_manual\_v1.12.pdf,
and open the example model "simple\_system.bngl". The manual describes v1.12;
features added since are described in this README and its release notes.

Source code is in the "src" directory. Example models are in the "models"
directory, with README files. Tests are in the "test" and "validate"
directories.

Enjoy your new network-free world!

## Building NFsim

### Linux and OSX

Make sure you have a recent version of CMake installed and run the following at
a terminal:

    mkdir build
    cd build
    cmake ..
    make

### Windows

Make sure you have a recent version of CMake, Cygwin (or MinGW), and Ninja
installed. Then run the following at a terminal:

    mkdir build
    cd build
    cmake -G "Ninja"
    ninja

### Optional: ExprTk Expression Parser

By default, NFsim uses muParser for mathematical expression parsing. You can
optionally build with ExprTk instead, which is a header-only, actively
maintained alternative. To enable ExprTk:

    mkdir build
    cd build
    cmake -DNFSIM_USE_EXPRTK=ON ..
    make

**Note:** ExprTk requires C++17, while the default muParser build uses C++11.

**Why ExprTk?**
- Header-only (no external dependencies)
- Actively maintained (muParser development ceased ~2016)
- Faster expression evaluation
- Future-proof against compiler evolution

### Optional: Build NFsim as a Library

By default, NFsim builds only the `NFsim` command-line executable. You can also
build a reusable `nfsim` library target from the same source tree:

    mkdir build
    cd build
    cmake -DNFSIM_BUILD_LIBRARY=ON ..
    make

To build only the library and skip the CLI executable:

    cmake -DNFSIM_BUILD_LIBRARY=ON -DNFSIM_BUILD_EXECUTABLE=OFF ..

To request a shared library instead of the default static library, also set:

    cmake -DNFSIM_BUILD_LIBRARY=ON -DBUILD_SHARED_LIBS=ON ..

The `NFSIM_USE_EXPRTK` option works with either the executable or library
targets.

If you've built it with Cygwin and you want to move or package up the
executable, you'll need to copy the the following DLLs along with it:

* cygwin1.dll
* cygstdc++-6.dll
* cygz.dll
* cyggcc_s-seh-1.dll.

## Practical Flags

- `-gml <integer>` sets the per-molecule-type maximum agent count.
- `-gml auto` (also accepts `none` and `nolimit`) disables this cap. This is
    useful for very large models where the default cap would stop the simulation.
- `-i <time>` sets the absolute simulation start time. The default is `0`; the
    value is used by `time()` and time-dependent functions and becomes the first
    value in the output time column. `-oTimes` values remain relative to the
    start of the requested simulation.
- A function that reads `time()` can be output with `-ogf`, but NFsim refuses
    a rate law that reads `time()`, directly or through another function. NFsim
    keeps each waiting time drawn at the last event, so such a rate would be
    re-read only when some reaction fires (RuleWorld/nfsim#112).

### Complex Bookkeeping and Molecularity

The command-line XML runner enforces the `+` in a rule on both sides by
default, as BioNetGen's network generator does:

- Two reactant patterns separated by `+` must match two different complexes.
    `A(b) + B(a) -> A(b!1).B(a!1)` does not close a ring inside one complex.
- A rule that deletes a bond and writes `+` between the two molecules it joined
    fires only if they end up in different complexes. `A(b!1).B(a!1) -> A(b) + B(a)`
    does not fire on a bond that leaves `A` and `B` connected some other way, and
    neither does `A(b!1).B(a!1) + C(s) -> A(b) + B(a) + C(s)`.
- A rule that deletes a bond and writes `.` between the two molecules it joined
    fires only if they are still in one complex afterwards. `A(b!1).B(a!1) -> A(b).B(a)`
    does not fire on a bond that is all that holds `A` and `B` together.
- The same two conditions hold for the molecules left behind when a rule deletes
    a molecule from the middle of a complex, as in
    `A(b!1).B(a!1,c!2).C(b!2) -> A(b) + C(b) DeleteMolecules`.
- A rule that writes molecules from different reactants in one product pattern
    without bonding them, such as `A(s~0) + B(a) -> A(s~1).B(a)`, never fires.
    NFsim prints a warning for it when the model is loaded.

Both checks need complex bookkeeping, so that is on by default as well.

- `-cb` and `-bscb` are still accepted and change nothing. `-cb` is the
    spelling emitted by BioNetGen's `complex=>1` setting.
- `-no-bscb` turns both checks off together. This was the behavior of a run
    given neither `-cb` nor `-bscb` before the checks became the default, and it
    is there for models written against that behavior. It cannot be combined
    with `-cb` or `-bscb`.
- With `-no-bscb`, a Species observable still auto-enables complex bookkeeping,
    because it is counted per complex.
- A program that calls `NFinput::initializeFromXML` directly chooses the policy
    through its `blockSameComplexBinding` argument, which has no default.

NFsim rejects `TotalRate` on Michaelis-Menten, Arrhenius, saturation, and local
function rate laws because those implementations do not define a compatible
macroscopic-rate interpretation.

## TFUN Support

NFsim accepts TFUN functions in the XML it consumes from BioNetGen-style
workflows. Current BioNetGen lowercase `tfun()` emits `<Function type="TFUN"
...>` XML, and NFsim supports that XML contract directly.

Supported TFUN XML forms:

- file-backed: `<Function type="TFUN" file="..." ctrName="..." method="linear|step">`
  The data file holds one point per line: the counter's value and the
  function's, separated by white space. Lines may end in LF, CRLF or CR. A
  `#` starts a comment that runs to the end of its line, so the header line
  of a BioNetGen `.tfun` file is read as one, and blank lines are skipped.
- inline: `<Function type="TFUN" mode="inline" ctrName="..." xData="..." yData="..." method="linear|step">`

Supported TFUN counters:

- `time` or `t`
- a model parameter
- an observable
- another global function

Supported TFUN expression placeholders:

- current BioNetGen-compatible placeholder: `__TFUN_VAL__`
- legacy placeholder retained for compatibility: `__TFUN__VAL__`

Compatibility notes:

- `linear` matches current BioNetGen lowercase `tfun()` output.
- `step` remains supported for older TFUN workflows.
- If the legacy placeholder `__TFUN__VAL__` is used without an explicit
  `method`, NFsim preserves the older default of `step`.
- TFUN values may appear inside larger expressions such as
  `(__TFUN_VAL__+5)/k_scale`.
- A rate law can read a table indexed by `time`. NFsim reads the table again
  at each of its x values, and where a `linear` table ramps, it places the next
  reaction where the integrated propensity reaches its random draw (#142). A
  rate law with arguments (a local function) that reads such a table is
  refused, because its value for each reactant is kept until that reactant
  changes.

This section describes the XML contract NFsim accepts. NFsim itself is not a
BNGL parser; BioNetGen remains responsible for translating `tfun()` or older
`TFUN()` syntax into XML or other NFsim inputs.

## Python Scripting Helper

For easier Python integration, a lightweight helper is provided at
`tools/python/nfsim_api.py`.

Example:

```python
from tools.python import run_nfsim, read_gdat

result = run_nfsim("path/to/model.xml", extra_args=["-sim", "100", "-oSteps", "100"])
rows = read_gdat(result.output_path)
print(rows[0])
```

## Optional Execution Profiler

NFsim can write a tab-separated CPU profile for investigating irregular reaction
work:

```text
NFsim -xml path/to/model.xml -sim 10 -oSteps 1 -profile profile.tsv
```

The report includes phase totals, per-reaction firing time, template matching,
membership updates, mapping operations, product preparation, connectivity
traversal by caller context, and bounded traversal-size distributions. Profiling
is opt-in and adds measurement overhead; ordinary runs do not write a profile.
Use fixed seeds and compare model output separately when evaluating a change.

## Release Notes

### v2.0.0 November, 2026

NFsim 2.0 is the first release since v1.14.3 and collects about fifty merged
pull requests. It adds compartmental and energy-based models, fixes a number of
rate errors, and changes some defaults, so results for existing models can
differ from v1.14.3. Numbers in parentheses are GitHub issues and pull requests.

**Changes that can alter results for existing models**

(a) `+` and `.` in rules are enforced by default, as in BioNetGen: two reactant
patterns joined by `+` must match different complexes, and a rule that breaks a
bond must leave the molecules separated (`+`) or still joined (`.`) as its
products say (#103, #51). `-cb` and `-bscb` are still accepted and change
nothing; `-no-bscb` restores the old behavior. Runs are slower because NFsim now
tracks complexes: up to about 30% on ordinary models and about 3x on TLBR.

(b) Several rates were wrong and are now corrected (#88, #89, #92): a rule with
a symmetric reactant pattern ran twice as fast as intended under every rate law
except the constant one; a TotalRate rule with a symmetric reaction center ran
too slowly (#90); a rule whose catalyst is a repeated structure ran too fast
(#87); a rule that breaks several bonds to open a ring could never fire; and the
Michaelis-Menten rate law lost precision when Km is small and the enzyme is in
excess (BioNetGen #323). A rule that deletes a whole species (no
`DeleteMolecules`) fired once per match of its pattern where BioNetGen gives it
one reaction per complex, so a homodimer went at twice the rule's rate and a
4-ring at four times; MatchOnce, which asks for the same counting, took no
effect on a reactant the rule changes, and was dropped with a warning on a
local-function rule (#148). A Michaelis-Menten rule cannot be corrected that
way, and stops with an error if it picks a complex that holds two matches of
such a reactant.

(c) Species observables report correct counts without a bookkeeping flag, and a
Species observable turns on complex tracking by itself (#48, #49, #89).

(d) All random draws in a simulation come from one generator owned by the
System (#67, #124, #125), so seeded runs give different trajectories from
v1.14.3. Builds of master between those changes picked reactants with the same
numbers that chose the rule, which biased which molecules reacted (#124).
Without `-seed`, NFsim now picks a seed, prints it, and uses it for the whole
run.

(e) NFsim keeps the total propensity as a running sum. After a fast rule
emptied, rounding could leave that sum at zero or at a multiple of a slow
rule's rate, so the slow rule never fired or fired at the wrong rate, and runs
using `-oTimes` never corrected it. NFsim also fired nothing while the total was
1e-9 or less, a rate in the model's own time unit. Where the same updates
repeat, rounding could also drift the sum steadily even though the total
never fell, starving a slow rule listed last. The sum is now recomputed once
it falls far below its recent peak and after every 4096 updates, and any
positive total can fire (#123, #130).

(f) The default molecule limit (`-gml`) is now 2^31 - 1 (#53), and the automatic
traversal limit is one more than the largest reactant pattern (#52).

**New features**

(g) Compartmental models (cBNGL). NFsim no longer stops on a compartments block.
It simulates molecules in volume and surface compartments, transport between
them (including MoveConnected), and volume scaling of rate constants (#43, #57,
#111). See the known issue on compartment sizes below.

(h) Energy-based models (eBNGL). Energy rules with Arrhenius rate laws are
expanded into ordinary rules when the model is loaded (#58, #100, #111).

(i) Fixed (clamped) species (#60), and the `include_products` and
`exclude_products` rule modifiers (#85).

(j) TFUN functions, file-backed or inline, with linear or step interpolation
(#69). See TFUN Support above.

(k) Functions can read `time()`, and `-i` sets the start time. A function that
reads `time()` can be output, but NFsim refuses a rate law that reads it,
directly or through another function, because such a rate was updated only when
a reaction fired (#112).

(l) For programs that embed NFsim: an optional library build
(`-DNFSIM_BUILD_LIBRARY=ON`), access to output observables, and a generator per
System so that several Systems can share a process (#67, #72, #73). An optional
ExprTk expression parser (`-DNFSIM_USE_EXPRTK=ON`) can replace muParser (#68).

(m) `-profile` writes a CPU profile of a run (#97).

(n) `-connect` is fast again on polymer models (#107, #114, #118), and several
of its bugs are fixed (#74, #120, #121).

(o) Performance work on reaction stepping, mapping sets, local functions and
energy-pattern rules (#44, #93, #94, #95, #98, #100).

**Disabled**

(p) RuleMonkey mode (`-rm`, `-rulemonkey`), added during this development
cycle, gave wrong results in some models and crashed on others (#116). NFsim
stops with an error if either flag is given.

**Other fixes**

(q) Out-of-bounds reads on a bond partner of the wrong molecule type (#113); a
crash on a negative composite rate law, now reported as an error (#119); rates
computed while a firing's products were outside the observables (#115, #117);
a use-after-free in `resetConcentrations()` (#81) and a double free at teardown
(#82); and many smaller fixes listed in the pull requests.

**Known issues**

- While complexes are tracked, which is the default and is required by Species
  observables, a reactant that a rule does not transform is counted once per
  complex. When a complex splits or joins
  farther from that reactant than the traversal limit reaches, the rule's
  propensity is not updated until the rule fires again or something else
  changes the reactant. A larger `-utl` does not fix this. For example, a rule
  `R(p!+) -> R(p!+) + Y()` acting on both ends of a chain that splits in the
  middle fires at about 70% of its correct rate. Check such models against a
  network-based simulation where one is feasible (#116).
- A compartment's size is fixed when the model is loaded; changing the
  parameter behind it later has no effect (#105).
- With `-connect`, a seeded run of some models is not identical to the same run
  without it, although averages agree (#116).

### v1.14.3 April, 2025

Changed compilation flags for MacOS to produce universal binary the runs on both arm64 and intel macs.

### v1.14.2 May, 2024

(a) Bugfix: Fixed a segmentation fault resulting from having a species start with zero concentration. Two simple models were added to the validation process to catch this error.
(b) Changes to the validation process were made. Previously, 15 trajectories would be generated for each model before testing the differences of NFsim and SSA versus NFsim and ODE. Now only one trajectory would be tested and more would generate upon faliure until a model failed 15 times in a row.

### v1.14.1 May, 2023

(a) Bugfix: Added a missing "\[" when no operations are present

### v1.14.0 February, 2023

Relabeling the version number to follow previous scheme to avoid further confusion. 

(a) Feature: New output format, `.nfevent.json`. JSON based format for all events that happen during a NFsim simulation. This format replaces the older reaction format and can be used with the `-rxnlog` command line argument. The use of `.nfevent.json` file extension is highly recommended since a JSON schema for the format is available. 
(b) Bugfix: Fixed an issue where `-connect` option would break with models that have remove operations.

### v1.2.2 April, 2022

(a) Bugfix: Disabled maxcputime option
(b) Bugfix: NFsim now checks for compartments and quits if it finds a compartments block
(c) Bugifx: Leftover debug statement was removed.

### v1.2.1 Feb, 2022

(a) Bugfix release. An accidentally leftover debug statement was removed.

### v1.2.0 Jan, 2022

(a) Added the original file-backed TFUN support, which can pull values from a file given a counter observable in a model.
(b) A new option is added to infer connectivity between reaction rules. This option allows the user to infer "connected" rules before a model is ran and each time a rule fires, only connected rules are checked for updates, instead of every possible rule in the system.

### v1.12.1 Aug, 2016

(a) Bugfix release. Addresses an error dealing with local function and species
labels. The error dealt with the way mappingSets where created and passed to
the localFunction evaluation. A test case (v19.bngl) was added that addresses
this case.

### v1.12   Dec, 2015

(a) Changes to how molecule instances are mapped to BasicRxn's
(BasicRxnClass::tryToAdd()). It was possible for certain kinds of rules that
the mappingSets were not updated correctly because the head molecule matched a
reactant pattern before and after a reaction event BUT the mapping was
different after the firing. In such cases, the mappingSet was not updated
properly. The logic was changed to fix this problem. (b) Further changes to
how molecule instances are mapped to a ReactionClass object (BasicRxn and
DORReactions). In particular it is often the case that graph symmetry leads to
a complex being able to map to a ReactionClass multiple times. Symmetry
considerations were being made on the reaction center but not on the context
components which led to an undercounting of the number of times a pattern agent
could match a rule instance in some edge cases (see v17.bngl in the validation
suite). This led to incorrect results or even NFSim crashes. NOTE: the changes
made in points (a) and (b) may cause some models to execute less efficiently.
(c) Fixed index bound checking in MoleculeType::getComponentStateName(). (d)
Updated the validation models to reflect current BNGL formatting standards.
Added a few new validation models that address the bugfixes included in this
version. This update also includes a Python version of the validation script.

### v1.11   Oct, 2012

(a) Molecules without components may be treated as population variables
rather than individual agents. This feature is useful for reducing memory
requirements when a simple molecule has very large population. A molecule
type may be flagged for treatment as a population variable by using the
"population" keyword following the molecule type definition in the BNGL
model file. See section 8.i in the documentation for further detail.
(b) Added a new Reaction Class called "FunctionProduct" that permits local
functions defined on two reactants. The local rate law must have the form
f(x)\*g(y), where x and y are tags on two distinct reactants. See section
7.c in the documentation for complete information.
(c) Fixed some problems evaluating complex-scoped local functions (CSLF). 
CSLF were not updated properly after reactions that split complexes or 
deleted molecules. As in v1.10, complex-scoped local functions are
enabled by default. A new command-line switch, -nocslf, has been added
which disables complex-scoped evaluation. 
(d) Improved efficiency for matching patterns with connected-to syntax
when the connected-to component does not have reaction center. This may
be especially notable in models with large complexes, e.g. polymerization.
		
### v1.10   Aug, 2011
        
(a) Command line parser now detects arguments that are not properly preceeded
by a dash, and generates a warning. (b) Includes a check when creating
template molecules that throws an error when users attempt to use Null or Trash
in reactant or observable patterns (anything that requires the creation of a
Template molecule). (c) fixed a bug introduced in v1.09 whereby a site was
allowed to bind to itself, for instance, in a dimerization rxn. (d) Support
for creating a new molecule bound to an existing molecule, as in a rule like
A(a) -> A(a!1).A(a!1). Existing code that implemented this feature did not
function properly with the check for null conditions before reactions were
fired. (e) Fixed bug in template molecule when clearing molecules after a
connected-to syntax search. In some cases, not all molecules were being
cleared, giving rise to situations where adding one observable created dangling
matches which affected the results of other observables. (f) NFsim is now
packaged with Network3, an updated version of the run\_network code to execute
ODE and SSA simulations. Network3 allows global functions in BNGL models among
other release features given here:
http://bionetgen.org/index.php/Release\_Notes. Note that Network3 does not
support On-The-Fly Stochastic Simulation (you will have to recompile Network2
to use this feature). (g)  Mac 32bit is no longer supported by NFsim, but you
can make executables for older Macs by recompiling the code on your own
machine. See the manual for instructions.

### v1.09   Apr, 2011

(a) NFsim now allows the mixing of integers and strings as component labels,
although if numbers and strings are mixed, all labels are parsed as strings,
NOT integers. Therefore, PLUS and MINUS keywords cannot be used if mixing
integer states and string states, and a warning will be generated if a state is
set to PLUS. One can only use PLUS or MINUS when ALL states are an integer
value greater than zero. This new behavior was needed to handle BNGL files that
used the convention of ~P specifying phosphorylated, and ~0 specifying
unphosphorylated.  (b) Fixed bug whereby if verbose option was turned on
without specifying an output file location, no output would be generated. Now,
output to a gdat file will be generated in these cases.  (c) Use of 'ss' input
argument to 'saveSpecies', which prints a a list of all species at the end of a
simulation, is now handled. This feature was implemented to allow future
support in BNG by restarting an NFsim simulation after it ends, which can be
done by parsing the output species list together with a BNGL model file. This
feature still has to be tested in BNG, and will likely be fully documented in
v1.10. The 'ss' flag writes the file to either system_name_nf.species, or a
file designated by the user.  (d) fixed memory leak in TemplateMolecules that
caused memory / performance issues with molecules having multiple identical
sites and a high degree of aggregation. (e) fixed csv error, where when the csv
flag is used, the header line is not comma delimited. (f)  nfsim now supports
intra-molecular binding. Previously these events were rejected as null events.

### v1.08 Dec, 2010

With the new TotalRate keyword, users are now able to specify whether or not to
use the microscopic (default) interpretation or macroscopic (TotalRate)
interpretation of rate laws. Now, NFsim convention matches BNG. Previously,
NFsim interpreted all rates as microscopic except for global functions, which
were interpreted as macroscopic. This is now also explained in the user manual.
Example models for the flagellar motor and oscillating gene expression have
been updated correspondingly so that they still produce the same results as in
the NFsim paper.

Users also now have the option of outputting gdat files in a comma delimited
format (csv), which makes parsing the output file easier in some circumstances,
using the flag "-csv". Additionally, a bug in the parameter scanning script was
fixed that caused the script to crash when scanning a model that includes the
local function syntax. 

### v1.07   Nov, 2010

A series of updates to the code were made in this release. (1) RNF files that
are not found produce an error message. Previously, no error message was given
and execution proceeded as if the RNF flag was not given. (2) Input flags to
NFsim can be given in the original format with a single dash (as in ./NFsim
-logo), or with the more "linuxy" style double dash (as in ./NFsim --logo).
(3) The parameter scan script had problems when parsing BNGL files with local
functions due to the '%' character. This is now fixed. (4) The universal
traversal limit is automatically set to be the size of the largest pattern in
the system, which can be overridden by passing the -utl flag. This allows
users who are unfamiliar with this speedup to still take advantage of it to a
certain extent. (5) the -rtag flag was added that allows NFsim to produce
output whenever a particular reaction, given by the -rtag flag, is given. This
allows, for instance, users to track the fates of single particles exactly
without using the comprehensive molecule output feature. (6) The above changes
are documented in an updated user manual.

### v1.06   Sept 28, 2010

Added scripts for running NFsim from Matlab, parameter scanning, and basic
parameter estimation. The manual is also updated to reflect these changes.
However, the precompiled executables of NFsim remain unchanged from v1.05, so
running them will give the old version number unless you recompile the code on
your own computer. Also, models that were used to compare the performance of
NFsim to DYNSTOC, RuleMonkey, and Kappa are now included with a readme file
under: models/performance\_test\_models.

### v1.052  

First publicly released stable build
