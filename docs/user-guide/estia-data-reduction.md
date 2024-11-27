# Polarization data reduction for ESTIA

Based on https://confluence.ess.eu/display/ESTIA/Polarised+Neutron+Reflectometry+%28PNR%29+-+Reduction+Notes

## Model

Intensity in the detector is related to the reflectivity of the sample by the model
```math
\begin{bmatrix}
I^{+} \\
I^{-}
\end{bmatrix}
\big(\lambda, j\big)
=
D(\lambda, j)
\begin{bmatrix}
1 - a^{\uparrow} & 1 - a^{\downarrow} \\
a^{\uparrow} & a^{\downarrow}
\end{bmatrix}
\begin{bmatrix}
(1 - f_2) & f_2 \\ f_2 & (1 - f_2)
\end{bmatrix}
\begin{bmatrix}
R^{\uparrow\uparrow} & R^{\downarrow\uparrow} \\
R^{\uparrow\downarrow} & R^{\downarrow\downarrow}
\end{bmatrix}
\begin{bmatrix}
(1 - f_1) & f_1 \\ f_1 & (1 - f_1)
\end{bmatrix}
\begin{bmatrix}
1 - p^{\uparrow} \\
1 - p^{\downarrow}
\end{bmatrix}
```
where

* $I^+$ is the intensity of the neutron beam transmitted by the analyzer
and $I^-$ is the intensity of the neutron beam reflected by the analyzer,
* $R^\cdot$ are the reflectivities of the sample,
  - $R^{\uparrow\uparrow}$ is the fraction of incoming neutrons with spin up that are reflected with spin up,
  - $R^{\uparrow\downarrow}$ is the fraction of incoming neutrons with spin up that are reflected with spin down,
  - etc..
* $a^\uparrow$ is the analyzer reflectivity for spin up neutrons and $a^\downarrow$ is the analyzer reflectivity for spin down neutrons,
* $p^\uparrow$ is the polarizer reflectivity for spin up neutrons and $p^\downarrow$ is the polarizer reflectivity for spin down neutrons,
* $f_1$ is the probability of spin flip by the polarizer spin flipper, $f_2$ is the probability of spin flip by the analyzer spin flipper
* $D$ represents the inhomogeneity from the beam- and detector efficiency (and all other polarization unrelated terms).

## Reducing a measurement

If the sample is measured at two different flipper settings $f_1=0, f_2=0$ and $f_1=1, f_2=0$, then we have four measurement in total:
```math
\begin{bmatrix}
I^{0+} \\
I^{0-}
\end{bmatrix}
\big(\lambda, j\big)
=
D(\lambda, j)
\begin{bmatrix}
1 - a^{\uparrow} & 1 - a^{\downarrow} \\
a^{\uparrow} & a^{\downarrow}
\end{bmatrix}
\begin{bmatrix}
R^{\uparrow\uparrow} & R^{\downarrow\uparrow} \\
R^{\uparrow\downarrow} & R^{\downarrow\downarrow}
\end{bmatrix}
\begin{bmatrix}
1 - p^{\uparrow} \\
1 - p^{\downarrow}
\end{bmatrix}
```
```math
\begin{bmatrix}
I^{1+} \\
I^{1-}
\end{bmatrix}
\big(\lambda, j\big)
=
D(\lambda, j)
\begin{bmatrix}
1 - a^{\uparrow} & 1 - a^{\downarrow} \\
a^{\uparrow} & a^{\downarrow}
\end{bmatrix}
\begin{bmatrix}
R^{\uparrow\uparrow} & R^{\downarrow\uparrow} \\
R^{\uparrow\downarrow} & R^{\downarrow\downarrow}
\end{bmatrix}
\begin{bmatrix}
1 - p^{\downarrow} \\
1 - p^{\uparrow}
\end{bmatrix}.
```

To simplify the above, collect the terms in the matrix $\mathbf{a}$
```math
\begin{bmatrix}
I^{0+} \\
I^{0-} \\
I^{1+} \\
I^{1-}
\end{bmatrix}
\big(\lambda, j\big)
=
D(\lambda, j)
\mathbf{a}(\lambda)
\begin{bmatrix}
R^{\uparrow\uparrow} \\
R^{\uparrow\downarrow} \\
R^{\downarrow\uparrow} \\
R^{\downarrow\downarrow}
\end{bmatrix}
(Q(\lambda, j)).
```

To compute the reflectivities, integrate over a region of (almost) constant $Q$
```math
\int_{Q\in[q_{n}, q_{n+1}]}
\mathbf{a}^{-1}(\lambda)
\begin{bmatrix}
I^{0+} \\
I^{0-} \\
I^{1+} \\
I^{1-}
\end{bmatrix}
\big(\lambda, j\big)
d\lambda dj
\approx
\int_{Q\in[q_{n}, q_{n+1}]}
D(\lambda, j)
d\lambda dj
\begin{bmatrix}
R^{\uparrow\uparrow} \\
R^{\uparrow\downarrow} \\
R^{\downarrow\uparrow} \\
R^{\downarrow\downarrow}
\end{bmatrix}
(q_{n+\frac{1}{2}}).
```
The integral on the righ-hand-side can be evaluated using the reference measurement.
$R$ was moved outside of the integral because if $Q$ is almost constant so is $R(Q)$.


### How to use the reference measurement to compute the integral over $D(\lambda, j)$?

For a reference measurement using flipper setting $f_1=0, f_2=0$ we have
```math
\begin{bmatrix}
I_{ref}^{+} \\
I_{ref}^{-}
\end{bmatrix}
\big(\lambda, j\big)
=
D(\lambda)
\begin{bmatrix}
1 - a^{\uparrow} & 1 - a^{\downarrow} \\
a^{\uparrow} & a^{\downarrow}
\end{bmatrix}
\begin{bmatrix}
R_{ref}^{\uparrow\uparrow} & R_{ref}^{\downarrow\uparrow} \\
R_{ref}^{\uparrow\downarrow} & R_{ref}^{\downarrow\downarrow}
\end{bmatrix}
\begin{bmatrix}
1 - p^{\uparrow} \\
1 - p^{\downarrow}
\end{bmatrix}
```
or, simplified by collecting the known terms in $r^{\pm}$
```math
\frac{1}{2}\bigg(\frac{I_{ref}^{+}(\lambda, j)}{r^+(\lambda, j)} + \frac{I_{ref}^{-}(\lambda, j)}{r^-(\lambda, j)}\bigg)
=
D(\lambda, j).
```
The expression for $D$ above can be used to evaluate integrals of $D$.

