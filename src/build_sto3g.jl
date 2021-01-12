"""
	build_sto3g(molecule::AbstractMolecule)

Build a STO-3G basis set for the given molecule.

This function only works for up to neon and uses data from EMSL.

This might change in the future, or even be deleted altogether.
In the future we will have basis set parsers, which will superseed this one.

Each function is given with ascending order of alpha value.
"""
function build_sto3g(molecule::AbstractMolecule)
	# Auxiliary function for checking which orbitals an element requires
	orbitalconfig(number) = number < 3 ? ["1s"] : ["1s", "2s", "2p"]

	T = typeof(molecule.coords[1])
	alphas_sto3g = Matrix{T}[]  # `alphas` for elements up to neon

	# H
	push!(alphas_sto3g, [0.3425250914e+01 0.6239137298e+00 0.1688554040e+00]')  # 1st shell
	# He
	push!(alphas_sto3g, [0.6362421394e+01 0.1158922999e+01 0.3136497915e+00]')
	# Li
	push!(alphas_sto3g, [0.1611957475e+02 0.2936200663e+01 0.7946504870e+00;    # 1st shell
						 0.6362897469e+00 0.1478600533e+00 0.4808867840e-01]')  # 2nd shell
	# Be
	push!(alphas_sto3g, [0.3016787069e+02 0.5495115306e+01 0.1487192653e+01;
						 0.1314833110e+01 0.3055389383e+00 0.9937074560e-01]')
	# B
	push!(alphas_sto3g, [0.4879111318e+02 0.8887362172e+01 0.2405267040e+01;
						 0.2236956142e+01 0.5198204999e+00 0.1690617600e+00]')
	# C
	push!(alphas_sto3g, [0.7161683735e+02 0.1304509632e+02 0.3530512160e+01;
						 0.2941249355e+01 0.6834830964e+00 0.2222899159e+00]')
	# N
	push!(alphas_sto3g, [0.9910616896e+02 0.1805231239e+02 0.4885660238e+01;
						 0.3780455879e+01 0.8784966449e+00 0.2857143744e+00]')
	# O
	push!(alphas_sto3g, [0.1307093214e+03 0.2380886605e+02 0.6443608313e+01;
						 0.5033151319e+01 0.1169596125e+01 0.3803889600e+00]')
	# F
	push!(alphas_sto3g, [0.1666791340e+03 0.3036081233e+02 0.8216820672e+01;
						 0.6464803249e+01 0.1502281245e+01 0.4885884864e+00]')
	# Ne
	push!(alphas_sto3g, [0.2070156070e+03 0.3770815124e+02 0.1020529731e+02;
						 0.8246315120e+01 0.1916266291e+01 0.6232292721e+00]')

	# Contraction coefficients that applies to all atoms
	coeffs_sto3g = [0.1543289673e+00 0.5353281423e+00 0.4446345422e+00;   # 1s
				   -0.9996722919e-01 0.3995128261e+00 0.7001154689e+00;   # 2s
					0.1559162750e+00 0.6076837186e+00 0.3919573931e+00]'  # all 2p

	sto3g = GaussianBasisFunction[]
	for i in 1:length(molecule)
		number = molecule.numbers[i]
		coord = molecule.coords[:, i]

		for orbital in orbitalconfig(number)
			if orbital == "1s"
				push!(sto3g, GaussianBasisFunction(
					coord,
					alphas_sto3g[number][:, 1],
					coeffs_sto3g[:, 1],
					0, 0, 0,
				))
			elseif orbital == "2s"
				push!(sto3g, GaussianBasisFunction(
					coord,
					alphas_sto3g[number][:, 2],
					coeffs_sto3g[:, 2],
					0, 0, 0,
				))
			elseif orbital == "2p"
				push!(sto3g, GaussianBasisFunction(
					coord,
					alphas_sto3g[number][:, 2],
					coeffs_sto3g[:, 3],
					1, 0, 0,
				))
				push!(sto3g, GaussianBasisFunction(
					coord,
					alphas_sto3g[number][:, 2],
					coeffs_sto3g[:, 3],
					0, 1, 0,
				))
				push!(sto3g, GaussianBasisFunction(
					coord,
					alphas_sto3g[number][:, 2],
					coeffs_sto3g[:, 3],
					0, 0, 1,
				))
			else
				error("Orbital $orbital not implemented for atomic number $number")
			end
		end
	end

	return sto3g
end
