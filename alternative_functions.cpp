// Calculate ratio R_SD
double calcRatio(const mat &detUp, const mat &detDown, double alpha, double beta, int i)
{
	double ratio = 0;
	double ri;
	if (i<nHalf)
	{
		for (int j = 0; j < nHalf; ++j)
		{
			ri = 0;
            for (j = 0; j < nDimensions; ++j) ri += r_new(i,j)*r_new(i,j);
            ri = sqrt(ri);
			ratio += phi(ri, alpha, j)*detUp[j][i];
		}
		return ratio;
	} else
	{
		for (int j = 0; j < nHalf; ++j)
		{
			ri = 0;
            for (j = 0; j < nDimensions; ++j) ri += r_new(i,j)*r_new(i,j);
            ri = sqrt(ri);
			ratio += phi(ri, alpha, j)*detDown[j][i-nHalf];
		}
		return ratio;
	}
}

// Update determinants
void update(mat &detUp, mat &detDown, int i, double ratio, double alpha, double beta)
{
	// Up determinant
	if (i<nHalf)
	{
		for (int j = 0; j < nHalf; ++j)
		{
			if (j!=i)
			{
				for (int k = 0; k < nHalf; ++k)
				{
					for (j = 0; j < nDimensions; ++j) ri += r_new(i,j)*r_new(i,j);
					ri = sqrt(ri);
					S_j[j]+=phi(ri, alpha, j)*detUp[k][j]
				}
			}
		}
		for (int j = 0; j < nHalf; ++j)
		{
			if (j!=i)
			{
				for (int k = 0; k < nHalf; ++k)
				{
					
					detUp[k][j]+=detUp[k][j]-S_j[j]*detUp[k][i]/ratio;
				}
			}
		}
		// i=j
		for (int k = 0; k < nHalf; ++k)
		{
			detUp[k][i]=detUp[k][i]/ratio;
		}
	}
	// Down determinant
	else
	{
		i = i - nHalf;
		for (int j = 0; j < nHalf; ++j)
		{
			if (j!=i)
			{
				for (int k = 0; k < nHalf; ++k)
				{
					for (j = 0; j < nDimensions; ++j) ri += r_new(i+nHalf,j)*r_new(i+nHalf,j);
					ri = sqrt(ri);
					S_j[j]+=phi(ri, alpha, j)*detDown[k][j]
				}
			}
		}
		for (int j = 0; j < nHalf; ++j)
		{
			if (j!=i)
			{
				for (int k = 0; k < nHalf; ++k)
				{
					
					detDown[k][j]+=detDown[k][j]-S_j[j]*detDown[k][i]/ratio;
				}
			}
		}
		// i=j
		for (int k = 0; k < nHalf; ++k)
		{
			detDown[k][i]=detDown[k][i]/ratio;
		}
		i = i + nHalf;
	}
}

// Jastrow factor ratio
double calculateJastrowRatio(mat &distance_old, mat &distance_new, double beta)
{
	double jastrowRatio;

	for (int i = 0; i < nParticles; ++i)
	{
		for (int j = 0; j < i; ++j)
		{
			jastrowRatio += distance_new[j][i]-distance_old[j][i];
		}
	}
	for (int i = 0; i < nParticles; ++i)
	{
		for (int j = i+1; j < nParticles; ++j)
		{
			jastrowRatio += distance_new[j][i]-distance_old[j][i];
		}
	}
	return jastrowRatio;
}

// Calculate distance
void calculateDistance(mat &distance, const mat r_old, double beta)
{
	for (int i = 0; i < nParticles; ++i)
	{
		for (int j = i+1; j < nParticles; ++j)
		{
			temp = diffR(r_old, i, j)
			//spin up
			if(((k < nHalf) && (l <nHalf)) || ((k>=nHalf && l>=nHalf)))
			{
				a=0.25;
				distance[k][l] = a*temp/(1+beta*temp);
			}
			//spin down
			else
			{
				a=0.5;
				distance[k][l] = a*temp/(1+beta*temp);
			}
		}
	}
}