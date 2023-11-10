# simulations

These files must be renamed and added into the `MC/CustomGenerators/PWGLF/` directory of [**AliDPG**](https://github.com/alisw/AliDPG):

* `GeneratorCustom.C`

These files must be added into the `EVGEN/` directory of [**AliRoot**](https://github.com/alisw/AliRoot):

* `AliGenSexaquarkReaction.cxx`
* `AliGenSexaquarkReaction.h`

Additionally, the class must also be added to the `EVGEN/CMakeLists.txt` and `EVGEN/EVGENLinkDef.h` files.

## Local Production

Once the files are added into [**AliRoot**](https://github.com/alisw/AliRoot) and [**AliDPG**](https://github.com/alisw/AliDPG), and these packages (besides [**AliPhysics**](https://github.com/alisw/AliPhysics)) are properly built and loaded, you need to execute this command to send a small production.

```
$ALIDPG_ROOT/bin/aliroot_dpgsim.sh --run 297595 --mode sim --uid 1 --nevents 1 --generator PWGLF:Hijing_Sexaquark:A1.8 --simulation SimulationDefaultIonTail --system Pb-Pb --detector Custom
```

## Testing

* **First successful attempt**

  1. Install latest tags of every package ([**AliRoot**](https://github.com/alisw/AliRoot), [**AliPhysics**](https://github.com/alisw/AliPhysics), [**AliDPG**](https://github.com/alisw/AliDPG), [**alidist**](https://github.com/alisw/alidist)) and add the generators as described above.
  2. Enter the environment:
     ```
     alienv enter AliRoot/latest,AliPhysics/latest,AliDPG/latest
     ```
  3. Download OCDB file via:
     ```
     alien.py cp alien:///alice/sim/2020/LHC20e3/OCDB/297595/OCDBsim.root file://OCDBsim.root
     ```
  4. Copy the current `DetectorCustom.C` macro into the test directory.
  5. Execute:
     ```
     $ALIDPG_ROOT/bin/aliroot_dpgsim.sh --run 297595 --mode sim --uid 1 --nevents 1 --generator PWGLF:Hijing_Sexaquark:A1.8 --simulation SimulationDefaultIonTail --system Pb-Pb --detector Custom
     ```

* **ALIROOT-8850**

  Description: *Pb-Pb, 5.02 TeV - Hijing + injected Hypernuclei anchored to pass3 of LHC18q and LHC18r, Geant3, minimum bias*
  1. Log into `lxplus.cern.ch` (not lxplus9!)
  2. Download OCDB file
     - load recent software (for example, `alienv enter O2Physics/daily-20231108-0100-1`)
     - download via `alien.py cp alien:///alice/sim/2020/LHC20e3/OCDB/297595/OCDBsim.root file://OCDBsim.root`
     - get out of env.
  3. Enter the environment:
     ```
     alienv enter AliPhysics/v5-09-54p-01_O2-1,AliDPG/prod-202305-01-1,jemalloc/latest
     ```
  4. Execute
     ```
     $ALIDPG_ROOT/bin/aliroot_dpgsim.sh --run 297595 --mode sim --uid 1 --nevents 1 --generator PWGLF:Hijing_Nuclex007 --simulation SimulationDefaultIonTail --system Pb-Pb
     ```
  Reference: any JDL file of the respective production tag LHC23f1.

