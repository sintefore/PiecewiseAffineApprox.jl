PiecewiseLinearApprox release notes
===================================

Version 0.5.0 (2024-05-31)
--------------------------
- Using affine rather than linear in naming
- Improve BigM calculation
- Pass arguments with Algorithm structs in favor of kwargs
- Add MIT license

Version 0.4.0 (2024-01-05)
--------------------------
- Rename package to PiecewiseAffineApprox to avoid name collision in General registry

Version 0.3.4 (2022-04-07)
--------------------------
- Fix type of pwlinear parameter for JuMP-variables

Version 0.3.3 (2022-03-10)
--------------------------
- Update dependencies to JuMP 0.23/1.0 (#13)
- Fix most broken tests

Version 0.3.2 (2022-01-28)
--------------------------
- Workaround for bigM with normals containing near-zeros (#10)
- Workaround for solver issue with strict=:none for compressor duty (#10)

Version 0.3.1 (2022-01-24)
--------------------------
- Improved calculation of big-M (#9)

Version 0.3.0 (2021-12-16)
--------------------------
- Add support for 2D functions (#3)
- New and breaking typed interface (#7)
- Add support for approximations based on Magnani/Boyd paper (#8)


Version 0.2.0 (2020-05-27)
--------------------------
- Added interpolation method (#1)
- Added optional plotting (#2)
- Added support for concave functions (#2)
- Breaking interface changess

Version 0.1.0 (2020-05-04)
--------------------------
- Initial release