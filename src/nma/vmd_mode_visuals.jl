"""
    write_vmd_mode(filename, coords, vecs, mode_idx, inv_M; n_frames=30)

Export a specific normal mode as a multi-frame XYZ trajectory for VMD/PyMOL.

This function creates an animation of a single vibrational mode by oscillating 
the structure along the eigenvector using a sine wave.

Parameters
----------
filename : String
    Path to the output .xyz file.
coords : Matrix{Float64}
    N×3 matrix of equilibrium Cartesian coordinates.
vecs : Matrix{Float64}
    The mass-weighted eigenvectors (mode shapes).
mode_idx : Int
    The index of the mode to visualize (e.g., 7 for the first vibrational mode).
inv_M : Diagonal{Float64}
    The inverse mass matrix (M⁻¹) or inv_sqrt_M (M⁻¹/²) to transform 
    the mode back to Cartesian space.
n_frames : Int, optional
    Number of frames in the animation (default: 30).

Notes
-----
- The transformation used is dₖ = M⁻¹/² * vₖ.
- The displacement is normalized to ensure the animation is visible.
- The oscillation follows the scaling factor: scale = A * sin(2π * frame / n_frames).
"""
function vmd_mode_visuals(filename, coords, vecs, mode_idx, inv_M; n_frames=30)
    N = size(coords, 1)

    # Transform the k-th eigenvector back to Cartesian displacement: dₖ = M⁻¹/² * vₖ
    mode = inv_M * vecs[:, mode_idx]
    
    # Normalize the displacement vector so the motion is consistent
    mode = mode / norm(mode)

    open(filename, "w") do io
        for frame in 1:n_frames
            # Create a smooth harmonic oscillation (sine wave): sin(theta)
            # A scale of 2.0 Å is usually a good starting point for visibility
            scale = 2.0 * sin(2π * frame / n_frames)

            # Standard XYZ format header
            println(io, N)
            println(io, "Mode $mode_idx - Frame $frame")  #println(io, "Mode $mode_idx")

            for i in 1:N
                # Extract x, y, z displacement components for atom i
                dx, dy, dz = mode[3i-2], mode[3i-1], mode[3i]
                
                # Write Atom Type and New Coordinates (Original + Oscillated)
                # Note: Using "C" as a generic atom type for XYZ readers
                println(io, "C $(coords[i,1] + scale*dx) $(coords[i,2] + scale*dy) $(coords[i,3] + scale*dz)")
            end
        end
    end
end

const write_vmd_mode = vmd_mode_visuals
