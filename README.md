# simulations

These files must be renamed and added into the `MC/CustomGenerators/PWGLF/` directory of [**AliDPG**](https://github.com/alisw/AliDPG):

* `GeneratorCustom.C`

These files must be added into the `EVGEN/` directory of [**AliRoot**](https://github.com/alisw/AliRoot):

* `AliGenSexaquarkReaction.cxx`
* `AliGenSexaquarkReaction.h`

Additionally, the class must also be added to the `EVGEN/CMakeLists.txt` and `EVGEN/EVGENLinkDef.h` files.

## Local Production

Once the files are added into [**AliRoot**](https://github.com/alisw/AliRoot) and [**AliDPG**](https://github.com/alisw/AliDPG), and these packages (besides [**AliPhysics**](https://github.com/alisw/AliPhysics)) are properly built and loaded, you need to execute this command to send a small production.

(**Note**: it's necessary to use the `OCDBsim.root` and `OCDBrec.root` callibration files that correspond to the run number you want, as described below in the **First successful attempt**)

```
$ALIDPG_ROOT/bin/aliroot_dpgsim.sh --run <run_number> --mode sim,rec,qa --uid 1 --nevents <n_events> --generator PWGLF:Hijing_Sexaquark:<reaction_channel><sexaquark_mass> --simulation SimulationDefaultIonTail --system Pb-Pb --detector NoAD
```
* `<run_number>` : run number
* `<n_events>` : number of events
* `<reaction_channel>` : reaction channel, could be
   - `A` : AntiSexaquark + Neutron -> AntiLambda, K0
   - `B` : AntiSexaquark + Neutron -> AntiLambda, K0, Pi-, Pi+
   - `C` : AntiSexaquark + Neutron -> AntiProton, K0, K0, Pi+
   - `D` : AntiSexaquark + Proton  -> AntiLambda, K+
   - `E` : AntiSexaquark + Proton  -> AntiLambda, K+, Pi-, Pi+
   - `F` : AntiSexaquark + Proton  -> AntiProton, K+, K0, Pi+
   - `G` : AntiSexaquark + Neutron -> Xi+, Pi-
   - `H` : AntiSexaquark + Proton  -> AntiProton, K+, K+, Pi0
* `<sexaquark_mass>` : hypothesized mass of the sexaquark, could go from 1.73 to 2.01 (in GeV/c^2)

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
  4. Execute:
     ```
     $ALIDPG_ROOT/bin/aliroot_dpgsim.sh --run 297595 --mode sim --uid 1 --nevents 1 --generator PWGLF:Hijing_Sexaquark:A1.8 --simulation SimulationDefaultIonTail --system Pb-Pb --detector NoAD
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

