//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
//  mctrueke                                                                    //
//  A multi-core implementation of the exchange Monte Carlo method.             //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
//  Copyright © 2015 Cristobal A. Navarro, Wei Huang.                           //
//                                                                              //
//  This file is part of mctrueke.                                              //
//  mctrueke is free software: you can redistribute it and/or modify            //
//  it under the terms of the GNU General Public License as published by        //
//  the Free Software Foundation, either version 3 of the License, or           //
//  (at your option) any later version.                                         //
//                                                                              //
//  mctrueke is distributed in the hope that it will be useful,                 //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of              //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               //
//  GNU General Public License for more details.                                //
//                                                                              //
//  You should have received a copy of the GNU General Public License           //
//  along with mctrueke.  If not, see <http://www.gnu.org/licenses/>.           //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////

#ifndef _CPUTOOLS_H_
#define _CPUTOOLS_H_

double computeF(const int *hlat, const int width, const int height, const int length){
		double F = 0.0;
		double mfx1=0.0;
		double mfx2=0.0;
		double mfy1=0.0;
		double mfy2=0.0;
		double mfz1=0.0;
		double mfz2=0.0;
		double k = 2.0*PI/(double)width;
		int idx, idy, idz;
		for(idx=0; idx<width; idx++){
			for(idy=0; idy<height; idy++){			
				for(idz=0; idz<length; idz++){					
					
						int i=idz*width*length+idy*width+idx;
						mfx1 += hlat[i] * cos(k * (double)idx);					
					    mfx2 += hlat[i] * sin(k * (double)idx);					
					
						mfy1 += hlat[i] * cos(k * (double)idy);					
					    mfy2 += hlat[i] * sin(k * (double)idy);					

						mfz1 += hlat[i] * cos(k * (double)idz);					
					    mfz2 += hlat[i] * sin(k * (double)idz);					

				}
			}
		}
		//printf("mfx1=%f mfx2=%f\nmfy1=%f mfy2=%f\nmfz1=%f mfz2=%f\n", mfx1, mfx2, mfy1, mfy2, mfz1, mfz2);
		//getchar();

		mfx1 *= mfx1;
		mfx2 *= mfx2;

		mfy1 *= mfy1;
		mfy2 *= mfy2;

		mfz1 *= mfz1;
		mfz2 *= mfz2;
		F = (1.0/(double)(3.0*width*length*height)) * ( mfx1 + mfx2 + mfy1 + mfy2 + mfz1 + mfz2 );
		return F;
}

void multicore_sweep(int *s, const int *H, const float h, const int L, const float B, const int tid, pcg32_random_t *rng, int sthreads){
    #pragma omp for
    for(int z = 0; z < L; z += 2){
        for(int y = 0; y < L; ++y){
            for(int x = 0; x < L; ++x){
                float ide = (float)s[C(x, y, z, L)] * (	((float)s[C((x==0)?L-1:x-1, y, z, L)] + s[C((x==L-1)?0:x+1, y, z, L)] + 
                                                                s[C(x, (y==0)?L-1:y-1, z, L)] + s[C(x, (y==L-1)?0:y+1, z, L)] +
                                                                s[C(x, y, (z==0)?L-1:z-1, L)] + s[C(x, y, (z==L-1)?0:z+1, L)])+ h*(float)H[C(x, y, z, L)] );

                //printf("x = %i   y = %i  z = %i\n", x, y, z); fflush(stdout);
                if(ide <= 0.0f || pcgrand(rng) < expf( ide * B ) ){
                //if(ide <= 0.0f || pcgrand(&(rng[tid])) < expf(ide * B) ){
                //if(ide <= 0.0f || 0.45 < expf( ide * B ) ){
                    s[C(x, y, z, L)] = -s[C(x, y, z, L)];
                }
            }
        }
    }
    #pragma omp barrier
    #pragma omp for
    for(int z = 1; z < L; z += 2){
        for(int y = 0; y < L; ++y){
            for(int x = 0; x < L; ++x){
                float ide = (float)s[C(x, y, z, L)] * (	((float)s[C((x==0)?L-1:x-1, y, z, L)] + s[C((x==L-1)?0:x+1, y, z, L)] + 
                                                                s[C(x, (y==0)?L-1:y-1, z, L)] + s[C(x, (y==L-1)?0:y+1, z, L)] +
                                                                s[C(x, y, (z==0)?L-1:z-1, L)] + s[C(x, y, (z==L-1)?0:z+1, L)])+ h*(float)H[C(x, y, z, L)] );

                //printf("x = %i   y = %i  z = %i\n", x, y, z); fflush(stdout);
                if(ide <= 0.0f || pcgrand(rng) < expf( ide * B ) ){
                //if(ide <= 0.0f || pcgrand(&(rng[tid])) < expf(ide * B) ){
                //if(ide <= 0.0f || 0.45 < expf( ide * B ) ){
                    s[C(x, y, z, L)] = -s[C(x, y, z, L)];
                }
            }
        }
    }
}

void cpu_sweep(int *s, const int *H, const float h, const int L, const float B, const int tid, pcg32_random_t *trng){
	for(int z = 0; z < L; ++z){
		for(int y = 0; y < L; ++y){
			for(int x = 0; x < L; ++x){
				float ide = (float)s[C(x, y, z, L)] * (	((float)s[C((x==0)?L-1:x-1, y, z, L)] + s[C((x==L-1)?0:x+1, y, z, L)] + 
																s[C(x, (y==0)?L-1:y-1, z, L)] + s[C(x, (y==L-1)?0:y+1, z, L)] +
																s[C(x, y, (z==0)?L-1:z-1, L)] + s[C(x, y, (z==L-1)?0:z+1, L)])+ h*(float)H[C(x, y, z, L)] );

				//printf("x = %i   y = %i  z = %i\n", x, y, z); fflush(stdout);
				if(ide <= 0.0f || pcgrand(trng) < expf( ide * B ) ){
				//if(ide <= 0.0f || pcgrand(&(rng[tid])) < expf(ide * B) ){
				//if(ide <= 0.0f || 0.45 < expf( ide * B ) ){
					s[C(x, y, z, L)] = -s[C(x, y, z, L)];
				}
			}
		}
	}
}

int computeM(const int *s, int N){
	int m = 0;
	for(int k = 0; k < N; ++k){
		m += s[k];
	}
	return m;
}

float computeE(const int *s, const int *H, float h, int L){
	float e = 0;
	for(int z = 0; z < L; ++z){
		for(int y = 0; y < L; ++y){
			for(int x = 0; x < L; ++x){
				e -= (float)s[C(x, y, z, L)] * ((float)(s[C((x==L-1)?0:x+1, y, z, L)] + s[C(x, (y==L-1)?0:y+1, z, L)] + s[C(x, y, (z==L-1)?0:z+1, L)])+ h*(float)H[C(x, y, z, L)]	);
			}
		}
	}
	return e;
}

// set a random Hi for the modle of random-field
void random_Hi(int *hat, int N, pcg32_random_t *trng){
for(int i=0; i < N; i++){
		if(pcgrand(trng) >= 0.5)
			hat[i] =  1;
		else
			hat[i] = -1;
	}	
}

#endif
