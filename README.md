# simulations

These files must be renamed and added into the `MC/CustomGenerators/PWGLF/` directory of [**AliDPG**](https://github.com/alisw/AliDPG):

* `GeneratorCustom.C`

These files must be added into the `EVGEN/` directory of [**AliRoot**](https://github.com/alisw/AliRoot):

* `AliGenSexaquarkReaction.cxx`
* `AliGenSexaquarkReaction.h`

Additionally, the class must be added to the `EVGEN/CMakeLists.txt` and `EVGEN/EVGENLinkDef.h` files.

## Local Production

Once the files are added into [**AliRoot**](https://github.com/alisw/AliRoot) and [**AliDPG**](https://github.com/alisw/AliDPG), and these packages (besides [**AliPhysics**](https://github.com/alisw/AliPhysics)) are properly built and loaded, you need to execute this command to send a small production.

```
$ALIDPG_ROOT/bin/aliroot_dpgsim.sh  --run 297595 --mode ocdb,full --uid 1 --generator PWGLF:Hijing_Sexaquark:A1.8 --nevents 1 --ocdb cvmfs --nocleanup
```
