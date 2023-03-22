# CT to MCNP

A script used to convert CT(nii) files to MCNP input files. Your must supply a tmol file which must contain `mode` card and `material` card. `examole.toml` is a template file.

# How to use

```
python main.py -c path/to/your/toml/file -d path/to/save/inp
```

# Toml Template

1. ct: path to your CT(nii) file(s).
2. Mode
   It is used for MCNP `mode` card. It is a string array.
3. Material
   It is used to convert CT images to voxels. You should type `[material.id]` for each material. `id` is the group and is used to generate `Mn` in MCNP.
   In each group, `hu_interval` indicates the CT value interval. `density` is the tissue density in this interval. `nucleon` and `fraction` is components and regarded fractions of this material. `nucleon` is a series of integers: Z * 1000 + A. `fraction` is a series of floating-point numbers. Their meanings is same as MCNP.
4. Source(Optinal)
   It is used to generate `SDEF` card. Use `keyword = [a series of values]`. If you want to use `SI` and `SP`, please use a dict shows in `example.toml`.
5. Tally(Optinal)
   It is used to generate `fmesh4` card. Only this tally card is supported by this script. `particle` indicates the type of particle recorded. And `FM `, `DE `, `DF` can also be added.
6. Output Control(Optinal)
   It is used to generate `nps`, `prdmp` and other card.
