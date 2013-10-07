/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <iomanip>
#include "graph.hpp"
#include "optics.hpp"
#include "source.hpp"
#include "io_timos.hpp"

#include <boost/timer/timer.hpp>

#include "mainloop.hpp"

template<class LoggerType,class RNG>int doOnePacket(const RunConfig& cfg,Packet pkt,
    LoggerType& logger,unsigned IDt,RNG& rng)
{
    unsigned Nstep=0;
    StepResult stepResult;
    Tetra currTetra = cfg.mesh.getTetra(IDt);
    Material currMat = cfg.mat[currTetra.matID];

	double dw,w0;
    float f_tmp[4] __attribute__((aligned(16)));
    float &n1 = f_tmp[0];
    float &n2 = f_tmp[1];
    float &ratio = f_tmp[2];
    unsigned IDt_next=IDt, IDm=currTetra.matID, IDm_next=IDm, IDm_bound;

    logger.eventLaunch(make_pair(pkt.p,pkt.d),IDt,1.0);

    // start another hop
    do
    {
        // draw a hop length; pkt.s = { physical distance, MFPs to go, time, 0 }
        ++Nstep;
        pkt.s = rng.draw_m128f1_log_u01();
        pkt.s = _mm_mul_ps(pkt.s,currMat.s_init);


        // attempt hop
        stepResult = currTetra.getIntersection(pkt.p,pkt.d,pkt.s);
        pkt.p      = stepResult.Pe;

        // loop while hitting a face in current step
        while (stepResult.idx != -1)
        {
            // extremely rarely, this can be a problem; we get no match in the getIntersection routine
            if(stepResult.idx > 3)
            {
                cerr << "Abnormal condition: stepResult.idx=" << stepResult.idx << ", IDte=" << stepResult.IDte << endl;
                return -1;
            }
            pkt.s = _mm_add_ps(pkt.s,_mm_mul_ps(stepResult.distance,currMat.s_prop));
            IDm_bound = cfg.mesh.getMaterial(stepResult.IDte);
            if (IDm == IDm_bound) { // no material change
                logger.eventBoundary(pkt.p,stepResult.IDfe,IDt,stepResult.IDte);
                IDt_next = stepResult.IDte;
            }
            else // boundary with material change
            {
                n2 = cfg.mat[IDm_bound].getn();
                n1 = currMat.getn();

                if (n1 == n2 || cfg.mat[IDm_bound].isMatched()) // no refractive index difference
                {
                    logger.eventBoundary(pkt.p,stepResult.IDfe,IDt,stepResult.IDte);
                    IDt_next = stepResult.IDte;
                }
                else {
                    logger.eventInterface(make_pair(pkt.p,pkt.d),stepResult.IDfe,stepResult.IDte);
                    __m128 Fn[4];

                    Fn[0] = currTetra.nx;
                    Fn[1] = currTetra.ny;
                    Fn[2] = currTetra.nz;
                    Fn[3] = _mm_setzero_ps();
                    _MM_TRANSPOSE4_PS(Fn[0],Fn[1],Fn[2],Fn[3]);
                    __m128 normal = Fn[stepResult.idx];
                    __m128 costheta = _mm_min_ps(_mm_set1_ps(1.0),_mm_sub_ps(_mm_setzero_ps(),dot3bcast(normal,pkt.d)));

                    ratio = n1/n2;
					__m128 n1_n2_ratio = _mm_load_ps(f_tmp);

                    __m128 sini_cosi_sint_cost = RefractSSE(n1_n2_ratio,costheta);

                    if (_mm_movemask_ps(_mm_cmplt_ss(_mm_set_ss(1.0),_mm_movehl_ps(sini_cosi_sint_cost,sini_cosi_sint_cost)))&1)
                    {
                        pkt.d = reflect(pkt.d,normal,sini_cosi_sint_cost);
                        logger.eventReflectInternal(pkt.p,pkt.d);
                    }
                    else {
                    __m128 d_p = _mm_add_ps(
									pkt.d,
									_mm_mul_ps(
										normal,
										costheta));

					    __m128 pr = FresnelSSE(n1_n2_ratio,sini_cosi_sint_cost);
                        if (_mm_movemask_ps(_mm_cmplt_ss(rng.draw_m128f1_u01(),pr))&1)
                        {
                            pkt.d = reflect(pkt.d,normal,sini_cosi_sint_cost);
                            logger.eventReflectFresnel(pkt.p,pkt.d);
                        }
                        else {
						    pkt.d = _mm_sub_ps(
							    _mm_mul_ps(
								    d_p,
								    _mm_shuffle_ps(n1_n2_ratio,n1_n2_ratio,_MM_SHUFFLE(2,2,2,2))),
							    _mm_mul_ps(
								    normal,
								    _mm_shuffle_ps(sini_cosi_sint_cost,sini_cosi_sint_cost,_MM_SHUFFLE(3,3,3,3))));
                            logger.eventRefract(pkt.p,pkt.d);
                            IDt_next = stepResult.IDte;
                        // configure material properties
                        } // if: fresnel reflection
                    }
                    pkt.a = getNormalTo(pkt.d);
                    pkt.b = cross(pkt.d,pkt.a);
                } // if: refractive index difference
            } // if: material change

            if (IDt_next == 0)
            {
                logger.eventExit(make_pair(pkt.p,pkt.d),stepResult.IDfe,pkt.w);
                return 0;
            }
            else if (IDt != IDt_next)
            {
                IDt = IDt_next;
                IDm_next = IDm_bound;
                currTetra = cfg.mesh.getTetra(IDt);
            }

            if (IDm != IDm_next)
            {
                IDm = IDm_next;
                currMat = cfg.mat[IDm];
                pkt.s = _mm_div_ss(_mm_movehdup_ps(pkt.s), _mm_set_ss(currMat.getMuT()));
            }
            stepResult=currTetra.getIntersection(pkt.p,pkt.d,pkt.s);
            pkt.p   = stepResult.Pe;
        }

        // stopped hitting faces: do drop/spin
        dw = currMat.getAbsorbedFraction()*pkt.w;
        logger.eventAbsorb(pkt.p,IDt,pkt.w,dw);
        pkt.w -= dw;

        // do roulette
        if (pkt.w < cfg.wmin && rng.draw_float_u01() < cfg.pr_win)
        {
            w0=pkt.w;
            pkt.w /= cfg.pr_win;
            logger.eventRouletteWin(w0,pkt.w);
        }

        // spin (new direction of travel)
        if (currMat.isScattering() && pkt.w >= cfg.wmin)
        {
			pkt = currMat.Scatter(pkt,rng.draw_float_u01(),rng.draw_m128f2_uvect());
            logger.eventScatter(pkt.d,pkt.d,currMat.getParam_g());
        }
    }
    while (pkt.w >= cfg.wmin && Nstep <= 10000);
    if(Nstep == 10000)
    {
        cerr << "Abnormal condition: packet retired after 10000 steps" << endl;
        cerr << "p=" << pkt.p << " d=" << pkt.d << " a=" << pkt.a << " b=" << pkt.b << endl;
        cerr << "IDt=" << IDt << " |d|=" << norm(pkt.d) << endl;
    }
    logger.eventDie(pkt.w);
    return 0;
}
