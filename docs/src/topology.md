# Topological Invariants

PhoXonic.jl provides functions to compute topological invariants of photonic and phononic band structures.

## Zak Phase (1D)

The Zak phase is the Berry phase acquired by a Bloch state when the wave vector k traverses the full Brillouin zone. For systems with inversion symmetry, the Zak phase is quantized to 0 or π.

```julia
lat = lattice_1d(1.0)
geo = Geometry(lat, Dielectric(1.0), [(Segment(0.2, 0.8), Dielectric(4.0))])
solver = Solver(Photonic1D(), geo, 128; cutoff=20)

result = compute_zak_phase(solver, 1:4; n_k=100)
result.phases  # Zak phase for each band
```

!!! note "Unit Cell Choice"
    Different choices of unit cell origin yield different Zak phases, but the difference between bands remains invariant.

## Wilson Loop Spectrum (2D)

The Wilson loop spectrum reveals band topology through the winding of eigenvalue phases. Non-zero winding indicates non-trivial topology.

```julia
lat = square_lattice(1.0)
geo = Geometry(lat, Dielectric(1.0), [(Circle([0.0, 0.0], 0.3), Dielectric(9.0))])
solver = Solver(TMWave(), geo, (32, 32); cutoff=5)

# Compute Wilson spectrum for bands 1-2
result = compute_wilson_spectrum(solver, 1:2; n_k_path=41, n_k_loop=50)

# Extract winding number
w = winding_number(result, 1)  # 0 = trivial, non-zero = topological
```

### Parameters

| Parameter | Description |
|-----------|-------------|
| `n_k_path` | Number of k-points along the scanning direction |
| `n_k_loop` | Number of k-points for Wilson loop integration |
| `loop_direction` | Direction of Wilson loop (`:b1` or `:b2`, default `:b2`) |

## API Reference

See [Advanced API - Topological Invariants](api-advanced.md#Topological-Invariants).

## References

- Zak, J. (1989). Berry's phase for energy bands in solids. *Phys. Rev. Lett.* 62, 2747.
- Yu, R. et al. (2011). Equivalent expression of Z2 topological invariant for band insulators using the non-Abelian Berry connection. *Phys. Rev. B* 84, 075119.
