tmp<fvVectorMatrix> tUEqn;

// Volumetric interpolation used, see references at top of RKSymFoam.C
const volVectorField gradpn
(
    fv::gaussGrad<scalar>::gradf(linear<scalar>(mesh).interpolate(pn), "gradpn")
    //fv::gaussGrad<scalar>::gradf(ghf*linear<scalar>(mesh).interpolate(rhok) + linear<scalar>(mesh).interpolate(pn), "gradpn")
);

/*//phig is required in buoyantBoussinesqPimpleFoam method
const volVectorField phig
(
    fv::gaussGrad<scalar>::gradf(ghf*linear<scalar>(mesh).interpolate(rhok), "phig")
);
*/

// Set up a UEqn matrix, containing:
// D: time derivative of U
// K: weighted sum of previous stage (explicit) convection and diffusion terms
// pnPredCoef weighted gradient of previous time step pressure 
// pnPredCoef = 0 equals Chorin method
// pnPredCoef = 1 equals Van Kan method
// possible fvOptions source term
// and finally, in case of implicit scheme:
// implicit convection and diffusion terms
if (scheme.explicitStage())
{
    tUEqn =
    (
        scheme.D(U)
      ==
        scheme.K()
      //-grad(p) + g*(rho/rho_ref)
      + scheme.C()*g*rhok
      - scheme.C()* gradpn*pnPredCoef
      + scheme.C()*fvOptions(U)
    );
}
else
{
    UEqnSourceA =
    (
      - fv::gaussConvectionScheme<vector>
        (
            mesh,
            phi,
            midPoint<vector>(mesh)
        ).fvmDiv(phi, U)
    );
	UEqnSourceB =
	    (
		fv::gaussLaplacianScheme<vector, scalar>
		(
		    mesh,
		    reverseLinear<scalar>(mesh),        // Scheme is not used
		    fv::uncorrectedSnGrad<vector>(mesh)
		).fvmLaplacian(nuf, U)
	      + fvc::surfaceIntegrate
		(
		    nuf
		  * reverseLinear<tensor>(mesh).dotInterpolate
		    (
			mesh.Sf(),
			dev2(Foam::T(
			    fv::gaussGrad<vector>::gradf(linear<vector>(mesh).interpolate(U), "gradU")
			))
		    )
		)
	    );

    tUEqn =
    (
        scheme.D(U)
      ==
        scheme.A()*UEqnSourceA
      + scheme.B()*UEqnSourceB
      + scheme.K()
      + scheme.C()*g*rhok
      - scheme.C()*gradpn*pnPredCoef
      + scheme.C()*fvOptions(U)
    );
}

fvVectorMatrix& UEqn = tUEqn.ref();

fvOptions.constrain(UEqn);
