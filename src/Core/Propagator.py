from .Matrix    import Matrix 
from .Precision import array_d
import math 
from numpy.linalg import inv, det

class Propagator : 
    def __init__(self, fieldmap= "", method ="") : 
        self.fieldmapfileName =""
        self.integrationMethod=""
    def Forward(self, interaction_point_index : int, win : Matrix , planes : list , tMat :Matrix = None  ) -> Matrix:
        if tMat !=None : 
            der = tMat.Inverse()    
        else :
            raise ValueError("Must pass propagation matrix input")  
        return Matrix( ( der.Transp() * (win * der )).matrix   )
    def Backward(self, interaction_point_index : int, win : Matrix , planes : list , tMat :Matrix = None  ) -> Matrix:
        if tMat !=None :
            der = tMat.Inverse()
        return Matrix( ( der.Transp() * (win * der )).matrix   )
    ### function doing stufff... have to go somewhere else.            
    def MultScat(self, scattering2 : float,  #no log term 
                       win   : Matrix, 
                       tx    : float= None, 
                       ty    : float = None  , 
                       dzMat : float = None , dzStep : float = None , backwardFit = False  ) -> Matrix  :
        if tx != None and ty != None :
            # is that already accounted for Outside the formula? 
            """
            const auto norm2cnoisetmp = norm2 * pow_2( std::min( std::abs( stv[4] ), iMeV ) );
            const auto radThick = t * std::sqrt( norm2 );
            // in a normal tracking run, for around 95% of cases, radThick is in the
            // interval [1e-7, 1e-1] - that's the region where the log turns nasty
            // because it decreases so rapidly, and it's not easily possible to
            // approximate it in a quick and dirty way
            //
            // be FMA friendly
            auto norm2cnoise = norm2cnoisetmp * radThick;   // sqrt(norm2) * norm2  * t          
            """
            norm2CNoise =  scattering2 * (1. + tx*tx + ty*ty )
            CoVTxTx =  ( 1. + tx*tx) * norm2CNoise
            CoVTyTy =  ( 1. + ty*ty) * norm2CNoise
            CoVTxTy =      tx *ty    * norm2CNoise
            C_noise_Mat =      array_d( [[ 0, 0, 0, 0, 0             ] ,
                                         [ 0, 0, 0, 0, 0             ] , 
                                         [ 0, 0 , CoVTxTx , CoVTxTy,0],
                                         [ 0, 0 , CoVTxTy,  CoVTyTy,0],
                                         [ 0, 0 ,    0   ,    0    ,0]])
            # Add extra terms for thick scatters 
            #see 
            #dzMat > 0 always, if backward propagation then it has to be -.5 factor
            UseThickness = True and dzMat != 0. 
            if UseThickness : 
                wallThickness    = dzMat 
                if backwardFit : factor = -.5
                else :           factor =   .5
                wallThicknessD   = wallThickness *  factor  #should be -.5 if MultScat bkwd
                thirds           = 1. / 3.
                wallThickness2_3 =  wallThickness *wallThickness * thirds            
                C_noise_Mat[0][0] = CoVTxTx * wallThickness2_3
                C_noise_Mat[1][0] = CoVTxTy * wallThickness2_3
                C_noise_Mat[0][1] = CoVTxTy * wallThickness2_3
                C_noise_Mat[1][1] = CoVTyTy * wallThickness2_3

                C_noise_Mat[2][0] = CoVTxTx * wallThicknessD
                C_noise_Mat[0][2] = CoVTxTx * wallThicknessD
                C_noise_Mat[2][1] = CoVTxTy * wallThicknessD
                C_noise_Mat[1][2] = CoVTxTy * wallThicknessD

                C_noise_Mat[3][0] = CoVTxTy * wallThicknessD
                C_noise_Mat[3][1] = CoVTyTy * wallThicknessD
                C_noise_Mat[0][3] = CoVTxTy * wallThicknessD
                C_noise_Mat[1][3] = CoVTyTy * wallThicknessD                
                #### maybe this ? 
                # state.covariance()( 0, 0 ) += 2 * dz * state.covariance()( 2, 0 ) + dz * dz * state.covariance()( 2, 2 );
                # state.covariance()( 2, 0 ) += dz * state.covariance()( 2, 2 );
                # state.covariance()( 1, 1 ) += 2 * dz * state.covariance()( 3, 1 ) + dz * dz * state.covariance()( 3, 3 );
                # state.covariance()( 3, 1 ) += dz * state.covariance()( 3, 3 );
            C_noise = Matrix(C_noise_Mat)
            """
                // now add the wall
                if ( applyScatteringCorrection ) {
                    m_scatteringTool->correctState( state, isept.material, thickness, upstream, pid );
                }
                if ( applyEnergyLossCorrection ) { dedxtool->correctState( state, isept.material, thickness, upstream, pid ); }

                // add the change in qOverP
                stateAtTarget.setQOverP( stateAtTarget.qOverP() + state.qOverP() - qop );

                // propagate the noise to the target. linear propagation, only
                // non-zero contributions
                const auto dz = ( upstream ? ztarget - z1 : ztarget - z2 );
                state.covariance()( 0, 0 ) += 2 * dz * state.covariance()( 2, 0 ) + dz * dz * state.covariance()( 2, 2 );
                state.covariance()( 2, 0 ) += dz * state.covariance()( 2, 2 );
                state.covariance()( 1, 1 ) += 2 * dz * state.covariance()( 3, 1 ) + dz * dz * state.covariance()( 3, 3 );
                state.covariance()( 3, 1 ) += dz * state.covariance()( 3, 3 );
                stateAtTarget.covariance() += state.covariance();
            """
            Idty = Matrix( array_d([[ 1. , 0 , 0, 0, 0 ],\
                                    [ 0 , 1. , 0, 0, 0 ],\
                                    [ 0 , 0 , 1., 0, 0 ],\
                                    [ 0 , 0 , 0, 1., 0 ],\
                                    [ 0 , 0 , 0, 0,  1.]]))
            #(1+WS)i.W
            # C_noise_Mat[0][0] +=  2*dzStep *C_noise_Mat[2][0] + dzStep*dzStep * C_noise_Mat[2][2]
            # C_noise_Mat[2][0] +=  dzStep   *C_noise_Mat[2][2] 
            # C_noise_Mat[0][2] +=  dzStep   *C_noise_Mat[2][2] 
            # C_noise_Mat[1][1] +=  2*dzStep *C_noise_Mat[3][1] + dzStep*dzStep * C_noise_Mat[3][3]
            # C_noise_Mat[3][1] +=  dzStep *C_noise_Mat[3][3]
            # C_noise_Mat[1][3] +=  dzStep *C_noise_Mat[3][3]
            to_invert = Idty + win*C_noise 
            if det( to_invert.matrix ) == 0.0 : 
                print("Determinant of inverting matrix is 0, will throw error")
            wout = ( Idty + win*C_noise ).Inverse() * win
        else : 
            raise ValueError("Invalid here")
            wout.AddMS_TX( scattering2)
            wout.AddMS_TY( scattering2)
        return wout    
        """
        I have the formulae for the mult. scatt. contribution for non-negligible tx and ty. 
        First, you have to rad length by sqrt(1+tx^2+ty^2) if the material plane is perp. to z (or, for any orientation, by 1/sin(a), 
        where a is the angle between the track and the plane). 
        Then you add to the error matrix the following terms:
        on (tx,tx) : (1+tx^2+ty^2)* (1+tx^2) * s2
        on (ty,ty) : (1+tx^2+ty^2)* (1+ty^2) * s2
        on (tx,ty) : (1+tx^2+ty^2)* tx * ty * s2
        where s2 = (14/p)^2.Lrad, where Lrad is the mumber of rad. lengths, accounting for the factor above.
        If you neglect the correlation terms, you van use the same functions as defore. Else you can code the matrix operation (1+WC)^(-1).W with the functions in Matrix5.cc  
        """


    # def PropagateTrackVector( before : np.array(), zbefore : float , zafter :float) :         
    #     """
    #     dx/dz = tx 
    #     dy/dz = ty 
    #     dtx/dz = [Q/P] * sqrt( 1+ tx2 + ty2 ) * (  ty* (txBx + Bz) - ( 1+tx2) *By )
    #     dty/dz = [Q/P] * sqrt( 1+ tx2 + ty2 ) * ( -tx* (tyBy + Bz) + ( 1+tx2) *By )
    #     """
    #     x   = before[0]
    #     y   = before[1]
    #     tx  = before[2]
    #     ty  = before[3]
    #     qop = before[4]


    # def BetheBlocCorrectionSimple( wIn , wOut ):
    #     """
    #         void StateSimpleBetheBlochEnergyCorrectionTool::correctState( LHCb::State& state, const Material* material,
    #                                                                 double wallThickness, bool upstream,
    #                                                                 LHCb::ParticleID ) const {
            
    #         double bbLoss = wallThickness * sqrt( 1. + std::pow( state.tx(), 2 ) + std::pow( state.ty(), 2 ) ) *
    #                             m_energyLossCorr * material->Z() * material->density() / material->A();
    #         bbLoss = std::min( m_maxEnergyLoss.value(), bbLoss );
    #         if ( !upstream ) bbLoss *= -1.;
    #         // apply correction - note for now only correct the state vector
    #         Gaudi::TrackVector& tX = state.stateVector();
    #         //  double minMomentumForEnergyCorrection = 10*Gaudi::Units::MeV;
    #         double qOverP        = 0.0;
    #         tX[4] < 0.0 ? qOverP = std::min( tX[4], -LHCb::Math::lowTolerance )
    #                                : qOverP = std::max( tX[4], LHCb::Math::lowTolerance );
    #         double newP         = 0.0;
    #         qOverP < 0.0 ? newP = std::min( 1.0 / qOverP - bbLoss, -m_minMomentumAfterEnergyCorr.value() )
    #                                          : newP = std::max( 1.0 / qOverP + bbLoss, m_minMomentumAfterEnergyCorr.value() );

    #         tX[4] = 1.0 / newP;

    #         // correction on cov
    #         if ( m_useEnergyLossError && m_sqrtEError > 0 ) {
    #         // error on dE is proportional to the sqrt of dE:
    #         // double sigmadE = m_sqrtEError * std::sqrt(std::abs(bbLoss))
    #         double                 err2 = m_sqrtEError * m_sqrtEError * std::abs( bbLoss ) * tX[4] * tX[4];
    #         Gaudi::TrackSymMatrix& cov  = state.covariance();
    #         cov( 4, 4 ) += err2;
    #         }
    #         }

    #     """
        