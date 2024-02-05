#include "GeoMagLib.h"

#include "EGM9615.h"

CDate::CDate()
{
	memset(this, 0x00, sizeof(CDate));
}

CDate::CDate(double sdate)
{
	memset(this, 0x00, sizeof(CDate));

	DecimalYear = sdate;
}

CDate::~CDate()
{

}

/* Sets WGS-84 parameters */
CEllipsoid::CEllipsoid()
{
    memset(this, 0x00, sizeof(CEllipsoid));

    a = 6378.137;				/*semi-major axis of the ellipsoid in */
    
	b = 6356.7523142;			/*semi-minor axis of the ellipsoid in */
    
	fla = 1 / 298.257223563;	/* flattening */
    
	eps = sqrt(1 - (b * b) / (a * a)); /*first eccentricity */
    epssq = (eps * eps);		/*first eccentricity squared */
    
	re = 6371.2;				/* Earth's radius */
}

CEllipsoid::~CEllipsoid()
{

}

CGeoid::CGeoid()
{
    memset(this, 0x00, sizeof(CGeoid));
}

/* Sets EGM-96 model file parameters */
CGeoid::CGeoid(const char p)
{
    NumbGeoidCols = 1441;	/* 360 degrees of longitude at 15 minute spacing */
    NumbGeoidRows = 721;	/* 180 degrees of latitude  at 15 minute spacing */
    NumbHeaderItems = 6;	/* min, max lat, min, max long, lat, long spacing*/

    ScaleFactor = 4;		/* 4 grid cells per degree at 15 minute spacing  */

    NumbGeoidElevs = NumbGeoidCols * NumbGeoidRows;

    Geoid_Initialized = 1;	/* TODO: Geoid will be initialized only if this is set to zero */

    /* If needed modify height referencing */
    switch (p)
    {
    
	case 'E':
    {
        break;
    }
	
	/* height above MSL */
    case 'M':
    {
        UseGeoid = 1;

        break;
    }
    
	default:
    {
        printf("\nProjection should be E or M");

        return;
    }
    }

}

CGeoid::~CGeoid()
{

}

double CGeoid::GetGeoidHeight(double lat, double lon)
{
    /* Latitude out of range */
    if ((lat < -90) || (lat > 90))
    {
        return -2;
    }

    OffsetY = (90.0 - lat) * ScaleFactor;

    /* Longitude out of range */
    if ((lon < -180) || (lon > 360))
    {
        return -1;
    }

    if (lon < 0.0)
    {
        OffsetX = (lon + 360.0) * ScaleFactor;
    }
    else
    {
        OffsetX = lon * ScaleFactor;
    }

    /*  Find Four Nearest Geoid Height Cells for specified Latitude, Longitude;   */
    /*  Assumes that (0,0) of Geoid Height Array is at Northwest corner:          */

    PostX = floor(OffsetX);

    if ((PostX + 1) == NumbGeoidCols)
    {
        PostX--;
    }

    PostY = floor(OffsetY);

    if ((PostY + 1) == NumbGeoidRows)
    {
        PostY--;
    }

    Index = (long)(PostY * NumbGeoidCols + PostX);

    ElevationNW = (double)GeoidHeightBuffer[Index];
    ElevationNE = (double)GeoidHeightBuffer[Index + 1];

    Index = (long)((PostY + 1) * NumbGeoidCols + PostX);

    ElevationSW = (double)GeoidHeightBuffer[Index];
    ElevationSE = (double)GeoidHeightBuffer[Index + 1];

    /*  Perform Bi-Linear Interpolation to compute Height above Ellipsoid:        */

    DeltaX = OffsetX - PostX;
    DeltaY = OffsetY - PostY;

    UpperY = ElevationNW + DeltaX * (ElevationNE - ElevationNW);
    LowerY = ElevationSW + DeltaX * (ElevationSE - ElevationSW);

    DeltaHeight = UpperY + DeltaY * (LowerY - UpperY);

    return DeltaHeight;
}

CCoordGeodetic::CCoordGeodetic()
{
    memset(this, 0x00, sizeof(CCoordGeodetic));
}

CCoordGeodetic::CCoordGeodetic(CGeoid* Geoid, double lat, double lon, double alt)
{
    memset(this, 0x00, sizeof(CCoordGeodetic));

    lambda = lon;
    
	phi = lat;
    
	HeightAboveGeoid = alt;

    /* This converts the height above mean sea level to height above the WGS-84 ellipsoid */
    if (Geoid->UseGeoid == 1)
    {
        double DeltaHeight = Geoid->GetGeoidHeight(lat, lon);

        HeightAboveEllipsoid = HeightAboveGeoid + DeltaHeight / 1000;
    }
    else
    {
        HeightAboveEllipsoid = HeightAboveGeoid;
    }
}

CCoordGeodetic::~CCoordGeodetic()
{

}

CCoordSpherical::CCoordSpherical()
{
    memset(this, 0x00, sizeof(CCoordSpherical));
}

CCoordSpherical::~CCoordSpherical()
{

}

void CCoordSpherical::GeodeticToSperical(CEllipsoid* Ellip, CCoordGeodetic* CoordGeodetic)
{
    /*
     ** Convert geodetic coordinates, (defined by the WGS-84
     ** reference ellipsoid), to Earth Centered Earth Fixed Cartesian
     ** coordinates, and then to spherical coordinates.
     */

    CosLat = cos(DEG2RAD(CoordGeodetic->phi));
    SinLat = sin(DEG2RAD(CoordGeodetic->phi));

    /* compute the local radius of curvature on the WGS-84 reference ellipsoid */

    rc = Ellip->a / sqrt(1.0 - Ellip->epssq * SinLat * SinLat);

    /* compute ECEF Cartesian coordinates of specified point (for longitude=0) */

    xp = (rc + CoordGeodetic->HeightAboveEllipsoid) * CosLat;
    zp = (rc * (1.0 - Ellip->epssq) + CoordGeodetic->HeightAboveEllipsoid) * SinLat;

    /* compute spherical radius and angle lambda and phi of specified point */

    r = sqrt(xp * xp + zp * zp);
    phig = RAD2DEG(asin(zp / r));		/* geocentric latitude */
    lambda = CoordGeodetic->lambda;		/* longitude */
}

CLegendreFunction::CLegendreFunction()
{
    memset(this, 0x00, sizeof(CLegendreFunction));
}

CLegendreFunction::CLegendreFunction(int n)
{
    memset(this, 0x00, sizeof(CLegendreFunction));

    NumTerms = n;

    Pcup = (double*)calloc((NumTerms + 1), sizeof(double));

    dPcup = (double*)calloc((NumTerms + 1), sizeof(double));
}

CLegendreFunction::~CLegendreFunction()
{
    if (Pcup)
    {
        free(Pcup);
    }

    if (dPcup)
    {
        free(dPcup);
    }
}

void CLegendreFunction::AssociatedLegendreFunction(CCoordSpherical* CoordSpherical, int nMax)
{
    double sin_phi = sin(DEG2RAD(CoordSpherical->phig)); /* sin  (geocentric latitude) */

    /* If nMax is less tha 16 or at the poles */
    if (nMax <= 16 || (1 - fabs(sin_phi)) < 1.0e-10)
    {
        CLegendreFunction::PcupLow(sin_phi, nMax);
    }
    else
    {
        CLegendreFunction::PcupHigh(sin_phi, nMax);
    }
}

void CLegendreFunction::PcupLow(double x, int nMax)
{
    int n, m, index, index1, index2, NumTerms;

    double k, z;

    double* schmidtQuasiNorm;

    Pcup[0] = 1.0;
    dPcup[0] = 0.0;

    /*sin (geocentric latitude) - sin_phi */
    z = sqrt((1.0 - x) * (1.0 + x));

    NumTerms = ((nMax + 1) * (nMax + 2) / 2);
    
	schmidtQuasiNorm = (double*)malloc((NumTerms + 1) * sizeof(double));

    /*	 First,	Compute the Gauss-normalized associated Legendre functions*/
    for (n = 1; n <= nMax; n++)
    {
        for (m = 0; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);
            
			if (n == m)
            {
                index1 = (n - 1) * n / 2 + m - 1;
                
				Pcup[index] = z * Pcup[index1];
                
				dPcup[index] = z * dPcup[index1] + x * Pcup[index1];
            }
            else if (n == 1 && m == 0)
            {
                index1 = (n - 1) * n / 2 + m;
                
				Pcup[index] = x * Pcup[index1];
                
				dPcup[index] = x * dPcup[index1] - z * Pcup[index1];
            }
            else if (n > 1 && n != m)
            {
                index1 = (n - 2) * (n - 1) / 2 + m;
                
				index2 = (n - 1) * n / 2 + m;
                
				if (m > n - 2)
                {
                    Pcup[index] = x * Pcup[index2];
                    
					dPcup[index] = x * dPcup[index2] - z * Pcup[index2];
                }
                else
                {
                    k = (double)(((n - 1) * (n - 1)) - (m * m)) / (double)((2 * n - 1) * (2 * n - 3));
                    
					Pcup[index] = x * Pcup[index2] - k * Pcup[index1];
                    
					dPcup[index] = x * dPcup[index2] - z * Pcup[index2] - k * dPcup[index1];
                }
            }
        }
    }
	
    /* Compute the ration between the the Schmidt quasi-normalized associated Legendre functions and the Gauss-normalized version. */
    schmidtQuasiNorm[0] = 1.0;
    
	for (n = 1; n <= nMax; n++)
    {
        index = (n * (n + 1) / 2);
        
		index1 = (n - 1) * n / 2;
        
		/* for m = 0 */
        schmidtQuasiNorm[index] = schmidtQuasiNorm[index1] * (double)(2 * n - 1) / (double)n;

        for (m = 1; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);
            
			index1 = (n * (n + 1) / 2 + m - 1);
            
			schmidtQuasiNorm[index] = schmidtQuasiNorm[index1] * sqrt((double)((n - m + 1) * (m == 1 ? 2 : 1)) / (double)(n + m));
        }
    }

    /* Converts the Gauss-normalized associated Legendre functions to
		the Schmidt quasi-normalized version using pre-computed relation
		stored in the variable schmidtQuasiNorm */

    for (n = 1; n <= nMax; n++)
    {
        for (m = 0; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);
            
			Pcup[index] = Pcup[index] * schmidtQuasiNorm[index];
            
            /* The sign is changed since the new WMM routines use derivative with respect to latitude insted of co-latitude */
			dPcup[index] = -dPcup[index] * schmidtQuasiNorm[index];
        }
    }

    if (schmidtQuasiNorm)
    {
        free(schmidtQuasiNorm);
    }
}

void CLegendreFunction::PcupHigh(double x, int nMax)
{
    double pm2, pm1, pmm, plm, rescalem, z, scalef;
    double* f1, * f2, * PreSqr;
    int k, kstart, m, n, NumTerms;

    NumTerms = ((nMax + 1) * (nMax + 2) / 2);

    if (fabs(x) == 1.0)
    {
        printf("Error in PcupHigh: derivative cannot be calculated at poles\n");
        //return FALSE;
        return;
    }

    f1 = (double*)malloc((NumTerms + 1) * sizeof(double));

    PreSqr = (double*)malloc((NumTerms + 1) * sizeof(double));

    f2 = (double*)malloc((NumTerms + 1) * sizeof(double));

    scalef = 1.0e-280;

    for (n = 0; n <= 2 * nMax + 1; ++n)
    {
        PreSqr[n] = sqrt((double)(n));
    }

    k = 2;

    for (n = 2; n <= nMax; n++)
    {
        k = k + 1;
        
		f1[k] = (double)(2 * n - 1) / (double)(n);
        
		f2[k] = (double)(n - 1) / (double)(n);
        
		for (m = 1; m <= n - 2; m++)
        {
            k = k + 1;
            
			f1[k] = (double)(2 * n - 1) / PreSqr[n + m] / PreSqr[n - m];
            
			f2[k] = PreSqr[n - m - 1] * PreSqr[n + m - 1] / PreSqr[n + m] / PreSqr[n - m];
        }
        
		k = k + 2;
    }

    /*z = sin (geocentric latitude) */
    z = sqrt((1.0 - x) * (1.0 + x));
    
	pm2 = 1.0;
    
	Pcup[0] = 1.0;
    
	dPcup[0] = 0.0;

    if (nMax == 0)
        //return FALSE;
        return;

    pm1 = x;
    
	Pcup[1] = pm1;
    
	dPcup[1] = z;
    
	k = 1;

    for (n = 2; n <= nMax; n++)
    {
        k = k + n;
        
		plm = f1[k] * x * pm1 - f2[k] * pm2;
        
		Pcup[k] = plm;
        
		dPcup[k] = (double)(n) * (pm1 - x * plm) / z;
        
		pm2 = pm1;
        pm1 = plm;
    }

    pmm = PreSqr[2] * scalef;
    
	rescalem = 1.0 / scalef;
    
	kstart = 0;

    for (m = 1; m <= nMax - 1; ++m)
    {
        rescalem = rescalem * z;

        /* Calculate Pcup(m, m)*/
        kstart = kstart + m + 1;
        
		pmm = pmm * PreSqr[2 * m + 1] / PreSqr[2 * m];
        
		Pcup[kstart] = pmm * rescalem / PreSqr[2 * m + 1];
        
		dPcup[kstart] = -((double)(m)*x * Pcup[kstart] / z);
        
		pm2 = pmm / PreSqr[2 * m + 1];
        
		/* Calculate Pcup(m+1, m)*/
        k = kstart + m + 1;
        
		pm1 = x * PreSqr[2 * m + 1] * pm2;
        
		Pcup[k] = pm1 * rescalem;
        
		dPcup[k] = ((pm2 * rescalem) * PreSqr[2 * m + 1] - x * (double)(m + 1) * Pcup[k]) / z;
        
		/* Calculate Pcup(n, m)*/
        for (n = m + 2; n <= nMax; ++n)
        {
            k = k + n;
            
			plm = x * f1[k] * pm1 - f2[k] * pm2;
            
			Pcup[k] = plm * rescalem;
            
			dPcup[k] = (PreSqr[n + m] * PreSqr[n - m] * (pm1 * rescalem) - (double)(n)*x * Pcup[k]) / z;
            
			pm2 = pm1;
            pm1 = plm;
        }
    }

    /* Calculate Pcup(nMax, nMax)*/
    rescalem = rescalem * z;
    
	kstart = kstart + m + 1;
    
	pmm = pmm / PreSqr[2 * nMax];

    Pcup[kstart] = pmm * rescalem;
    
	dPcup[kstart] = -(double)(nMax)*x * Pcup[kstart] / z;

    free(f1);
    free(PreSqr);
    free(f2);
}

CMagneticModel::CMagneticModel()
{
    memset(this, 0x00, sizeof(CMagneticModel));
}

CMagneticModel::CMagneticModel(int n)
{
    memset(this, 0x00, sizeof(CMagneticModel));

    num_terms = n;

    Main_Field_Coeff_G = (double*)calloc((num_terms + 1), sizeof(double));
    Main_Field_Coeff_H = (double*)calloc((num_terms + 1), sizeof(double));

    Secular_Var_Coeff_G = (double*)calloc((num_terms + 1), sizeof(double));
    Secular_Var_Coeff_H = (double*)calloc((num_terms + 1), sizeof(double));
}

CMagneticModel::CMagneticModel(const char* f)
{
    memset(this, 0x00, sizeof(CMagneticModel));

    if (f == nullptr)
    {
        printf("\nMagnetic Model file name not provided");

        return;
    }

    err = fopen_s(&infile, f, "rt");

    if (err)
    {
        printf("\nMagnetic Model file not found:%s", f);

        return;
    }

    // skipping the header while checking for end-of-file
    if (fgets(line, MAXLINELENGTH, infile) == eof)
    {
        printf("\nMagnetic Model file is empty");

        return;
    }

    eof = fgets(line, MAXLINELENGTH, infile);

    while (eof != nullptr)
    {
        // checking first value
        a = sscanf_s(line, "%d", &n);

        if (n > 12)
        {
            break;
        }

        if (n > nMax)
        {
            nMax = n;
        }

        eof = fgets(line, MAXLINELENGTH, infile);
    }

    fclose(infile);

    nMaxSecVar = nMax;

    //num_terms = CALCULATE_NUMTERMS(nMax);
    num_terms = (nMax * (nMax + 1) / 2 + nMax);

    Main_Field_Coeff_G = (double*)calloc((num_terms + 1), sizeof(double));
    Main_Field_Coeff_H = (double*)calloc((num_terms + 1), sizeof(double));

    Secular_Var_Coeff_G = (double*)calloc((num_terms + 1), sizeof(double));
    Secular_Var_Coeff_H = (double*)calloc((num_terms + 1), sizeof(double));

    CoefficientFileEndDate = epoch + 5;



    err = fopen_s(&infile, f, "rt");

    Main_Field_Coeff_H[0] = 0.0;
    Main_Field_Coeff_G[0] = 0.0;
    
	Secular_Var_Coeff_H[0] = 0.0;
    Secular_Var_Coeff_G[0] = 0.0;

    eof = fgets(line, MAXLINELENGTH, infile);

    sscanf(line, "%lf%s", &epoch, ModelName);
    //MagneticModel->epoch = epoch;

    eof = fgets(line, MAXLINELENGTH, infile);

    while (eof != nullptr)
    {
        if (line[0] == '9')
        {
            break;
        }

        sscanf(line, "%d%d%lf%lf%lf%lf", &n, &m, &gnm, &hnm, &dgnm, &dhnm);

        if (m <= n)
        {
            index = (n * (n + 1) / 2 + m);

            Main_Field_Coeff_G[index] = gnm;
            
			Secular_Var_Coeff_G[index] = dgnm;
            
			Main_Field_Coeff_H[index] = hnm;
            
			Secular_Var_Coeff_H[index] = dhnm;
        }

        eof = fgets(line, MAXLINELENGTH, infile);
    }

    CoefficientFileEndDate = epoch + 5;

    fclose(infile);
}

CMagneticModel::~CMagneticModel()
{
    if (Main_Field_Coeff_G)
    {
        free(Main_Field_Coeff_G);
    }

    if (Main_Field_Coeff_H)
    {
        free(Main_Field_Coeff_H);
    }

    if (Secular_Var_Coeff_G)
    {
        free(Secular_Var_Coeff_G);
    }

    if (Secular_Var_Coeff_H)
    {
        free(Secular_Var_Coeff_H);
    }
}

void CMagneticModel::TimelyModifyMagneticModel(CDate* UserDate, CMagneticModel* MagneticModel)
{
    int n, m, index, a, b;

    EditionDate = MagneticModel->EditionDate;

    epoch = MagneticModel->epoch;

    nMax = MagneticModel->nMax;

    nMaxSecVar = MagneticModel->nMaxSecVar;

    a = nMaxSecVar;

    b = (a * (a + 1) / 2 + a);

    strcpy_s(ModelName, 32, MagneticModel->ModelName);

    for (n = 1; n <= MagneticModel->nMax; n++)
    {
        for (m = 0; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);

            if (index <= b)
            {
                Main_Field_Coeff_H[index] = MagneticModel->Main_Field_Coeff_H[index] + (UserDate->DecimalYear - MagneticModel->epoch) * MagneticModel->Secular_Var_Coeff_H[index];
                Main_Field_Coeff_G[index] = MagneticModel->Main_Field_Coeff_G[index] + (UserDate->DecimalYear - MagneticModel->epoch) * MagneticModel->Secular_Var_Coeff_G[index];

                /* We need a copy of the secular var coef to calculate secular change */
                Secular_Var_Coeff_H[index] = MagneticModel->Secular_Var_Coeff_H[index];
                Secular_Var_Coeff_G[index] = MagneticModel->Secular_Var_Coeff_G[index];
            }
            else
            {
                Main_Field_Coeff_H[index] = MagneticModel->Main_Field_Coeff_H[index];
                Main_Field_Coeff_G[index] = MagneticModel->Main_Field_Coeff_G[index];
            }
        }
    }
}

CSphericalHarmonicVariables::CSphericalHarmonicVariables()
{
    memset(this, 0x00, sizeof(CSphericalHarmonicVariables));
}

CSphericalHarmonicVariables::CSphericalHarmonicVariables(int n)
{
    memset(this, 0x00, sizeof(CSphericalHarmonicVariables));

    nMax = n;

    RelativeRadiusPower = (double*)calloc((nMax + 1), sizeof(double));

    cos_mlambda = (double*)calloc((nMax + 1), sizeof(double));

    sin_mlambda = (double*)calloc((nMax + 1), sizeof(double));
}

CSphericalHarmonicVariables::~CSphericalHarmonicVariables()
{
    if (RelativeRadiusPower)
    {
        free(RelativeRadiusPower);
    }

    if (cos_mlambda)
    {
        free(cos_mlambda);
    }

    if (sin_mlambda)
    {
        free(sin_mlambda);
    }
}

void CSphericalHarmonicVariables::ComputeSphericalHarmonicVariables(CEllipsoid* Ellip, CCoordSpherical* CoordSpherical, int nMax)
{
    double cos_lambda, sin_lambda;
    int m, n;

    cos_lambda = cos(DEG2RAD(CoordSpherical->lambda));
    sin_lambda = sin(DEG2RAD(CoordSpherical->lambda));

    /* for n = 0 ... model_order, compute (Radius of Earth / Spherical radius r)^(n+2)
    for n  1..nMax-1 (this is much faster than calling pow MAX_N+1 times).      */
    RelativeRadiusPower[0] = (Ellip->re / CoordSpherical->r) * (Ellip->re / CoordSpherical->r);

    for (n = 1; n <= nMax; n++)
    {
        RelativeRadiusPower[n] = RelativeRadiusPower[n - 1] * (Ellip->re / CoordSpherical->r);
    }

    /*
     Compute cos(m*lambda), sin(m*lambda) for m = 0 ... nMax
           cos(a + b) = cos(a)*cos(b) - sin(a)*sin(b)
           sin(a + b) = cos(a)*sin(b) + sin(a)*cos(b)
     */
    cos_mlambda[0] = 1.0;
    sin_mlambda[0] = 0.0;

    cos_mlambda[1] = cos_lambda;
    sin_mlambda[1] = sin_lambda;

    for (m = 2; m <= nMax; m++)
    {
        cos_mlambda[m] = cos_mlambda[m - 1] * cos_lambda - sin_mlambda[m - 1] * sin_lambda;
        sin_mlambda[m] = cos_mlambda[m - 1] * sin_lambda + sin_mlambda[m - 1] * cos_lambda;
    }
}

CMagneticResults::CMagneticResults()
{
    memset(this, 0x00, sizeof(CMagneticResults));
}

CMagneticResults::~CMagneticResults()
{

}

void CMagneticResults::Summation(CLegendreFunction* LegendreFunction, CMagneticModel* MagneticModel, CSphericalHarmonicVariables* SphVariables, CCoordSpherical* CoordSpherical)
{
    int m, n, index;
    double cos_phi;

    Bz = 0.0;
    By = 0.0;
    Bx = 0.0;
	
    for (n = 1; n <= MagneticModel->nMax; n++)
    {
        for (m = 0; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);

            /*		    nMax  	(n+2) 	  n     m            m           m
                    Bz =   -SUM (a/r)   (n+1) SUM  [g cos(m p) + h sin(m p)] P (sin(phi))
                                    n=1      	      m=0   n            n           n  */
            
			/* Equation 12 in the WMM Technical report.  Derivative with respect to radius.*/
            Bz -= SphVariables->RelativeRadiusPower[n] *
                (MagneticModel->Main_Field_Coeff_G[index] * SphVariables->cos_mlambda[m] +
                    MagneticModel->Main_Field_Coeff_H[index] * SphVariables->sin_mlambda[m])
                * (double)(n + 1) * LegendreFunction->Pcup[index];

            /*		  1 nMax  (n+2)    n     m            m           m
                    By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                               n=1             m=0   n            n           n  */
            
			/* Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius. */
            By += SphVariables->RelativeRadiusPower[n] *
                (MagneticModel->Main_Field_Coeff_G[index] * SphVariables->sin_mlambda[m] -
                    MagneticModel->Main_Field_Coeff_H[index] * SphVariables->cos_mlambda[m])
                * (double)(m)*LegendreFunction->Pcup[index];
            
			/*		   nMax  (n+2) n     m            m           m
                    Bx = - SUM (a/r)   SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                               n=1         m=0   n            n           n  */
            
			/* Equation 10  in the WMM Technical report. Derivative with respect to latitude, divided by radius. */
            Bx -= SphVariables->RelativeRadiusPower[n] *
                (MagneticModel->Main_Field_Coeff_G[index] * SphVariables->cos_mlambda[m] +
                    MagneticModel->Main_Field_Coeff_H[index] * SphVariables->sin_mlambda[m])
                * LegendreFunction->dPcup[index];
        }
    }

    cos_phi = cos(DEG2RAD(CoordSpherical->phig));
    
	/* Special calculation for component - By - at Geographic poles.
		If the user wants to avoid using this function,  please make sure that
		the latitude is not exactly +/-90. An option is to make use the function
		MAG_CheckGeographicPoles.
	*/
	if (fabs(cos_phi) > 1.0e-10)
    {
        By = By / cos_phi;
    }
    else
    {
        CMagneticResults::SummationSpecial(MagneticModel, SphVariables, CoordSpherical);
    }
}

void CMagneticResults::SummationSpecial(CMagneticModel* MagneticModel, CSphericalHarmonicVariables* SphVariables, CCoordSpherical* CoordSpherical)
{
    int n, index;
    double k, sin_phi, * PcupS, schmidtQuasiNorm1, schmidtQuasiNorm2, schmidtQuasiNorm3;

    PcupS = (double*)malloc((MagneticModel->nMax + 1) * sizeof(double));

    PcupS[0] = 1;
    
	schmidtQuasiNorm1 = 1.0;

    By = 0.0;
    
	sin_phi = sin(DEG2RAD(CoordSpherical->phig));

    for (n = 1; n <= MagneticModel->nMax; n++)
    {
		/*Compute the ration between the Gauss-normalized associated Legendre
			functions and the Schmidt quasi-normalized version. This is equivalent to
			sqrt((m==0?1:2)*(n-m)!/(n+m!))*(2n-1)!!/(n-m)!
		*/
        index = (n * (n + 1) / 2 + 1);
        
		schmidtQuasiNorm2 = schmidtQuasiNorm1 * (double)(2 * n - 1) / (double)n;
        schmidtQuasiNorm3 = schmidtQuasiNorm2 * sqrt((double)(n * 2) / (double)(n + 1));
        schmidtQuasiNorm1 = schmidtQuasiNorm2;
        
		if (n == 1)
        {
            PcupS[n] = PcupS[n - 1];
        }
        else
        {
            k = (double)(((n - 1) * (n - 1)) - 1) / (double)((2 * n - 1) * (2 * n - 3));
            
			PcupS[n] = sin_phi * PcupS[n - 1] - k * PcupS[n - 2];
        }

        /*		  1 nMax  (n+2)    n     m            m           m
                By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                           n=1             m=0   n            n           n  */
        
		/* Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius. */
        By += SphVariables->RelativeRadiusPower[n] *
            (MagneticModel->Main_Field_Coeff_G[index] * SphVariables->sin_mlambda[1] -
                MagneticModel->Main_Field_Coeff_H[index] * SphVariables->cos_mlambda[1])
            * PcupS[n] * schmidtQuasiNorm3;
    }

    if (PcupS)
    {
        free(PcupS);
    }
}

/*This Function sums the secular variation coefficients to get the secular variation of the Magnetic vector.
INPUT :  LegendreFunction
		 MagneticModel
		 SphVariables
		 CoordSpherical
OUTPUT : MagneticResults

CALLS : MAG_SecVarSummationSpecial
 */
void CMagneticResults::SecVarSummation(CLegendreFunction* LegendreFunction, CMagneticModel* MagneticModel, CSphericalHarmonicVariables* SphVariables, CCoordSpherical* CoordSpherical)
{
    int m, n, index;
    double cos_phi;

    //MagneticModel->SecularVariationUsed = TRUE;
    MagneticModel->SecularVariationUsed = 1;

    Bz = 0.0;
    By = 0.0;
    Bx = 0.0;

    for (n = 1; n <= MagneticModel->nMaxSecVar; n++)
    {
        for (m = 0; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);

            /*		    nMax  	(n+2) 	  n     m            m           m
                    Bz =   -SUM (a/r)   (n+1) SUM  [g cos(m p) + h sin(m p)] P (sin(phi))
                                    n=1      	      m=0   n            n           n  */
            
			/*  Derivative with respect to radius.*/
            Bz -= SphVariables->RelativeRadiusPower[n] *
                (MagneticModel->Secular_Var_Coeff_G[index] * SphVariables->cos_mlambda[m] +
                    MagneticModel->Secular_Var_Coeff_H[index] * SphVariables->sin_mlambda[m])
                * (double)(n + 1) * LegendreFunction->Pcup[index];

            /*		  1 nMax  (n+2)    n     m            m           m
                    By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                               n=1             m=0   n            n           n  */
            
			/* Derivative with respect to longitude, divided by radius. */
            By += SphVariables->RelativeRadiusPower[n] *
                (MagneticModel->Secular_Var_Coeff_G[index] * SphVariables->sin_mlambda[m] -
                    MagneticModel->Secular_Var_Coeff_H[index] * SphVariables->cos_mlambda[m])
                * (double)(m)*LegendreFunction->Pcup[index];
				
            /*		   nMax  (n+2) n     m            m           m
                    Bx = - SUM (a/r)   SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                               n=1         m=0   n            n           n  */
                               
			/* Derivative with respect to latitude, divided by radius. */
            Bx -= SphVariables->RelativeRadiusPower[n] *
                (MagneticModel->Secular_Var_Coeff_G[index] * SphVariables->cos_mlambda[m] +
                    MagneticModel->Secular_Var_Coeff_H[index] * SphVariables->sin_mlambda[m])
                * LegendreFunction->dPcup[index];
        }
    }
	
    cos_phi = cos(DEG2RAD(CoordSpherical->phig));
    
	/* Special calculation for component By at Geographic poles */
	if (fabs(cos_phi) > 1.0e-10)
    {
        By = By / cos_phi;
    }
    else
    {
        CMagneticResults::SecVarSummationSpecial(MagneticModel, SphVariables, CoordSpherical);
    }
}

/*Special calculation for the secular variation summation at the poles.

INPUT:  MagneticModel
		SphVariables
		CoordSpherical
OUTPUT: MagneticResults
CALLS : none
 */
void CMagneticResults::SecVarSummationSpecial(CMagneticModel* MagneticModel, CSphericalHarmonicVariables* SphVariables, CCoordSpherical* CoordSpherical)
{
    int n, index;
    double k, sin_phi, * PcupS, schmidtQuasiNorm1, schmidtQuasiNorm2, schmidtQuasiNorm3;

    PcupS = (double*)malloc((MagneticModel->nMaxSecVar + 1) * sizeof(double));

    PcupS[0] = 1;

    schmidtQuasiNorm1 = 1.0;

    By = 0.0;

    sin_phi = sin(DEG2RAD(CoordSpherical->phig));

    for (n = 1; n <= MagneticModel->nMaxSecVar; n++)
    {
        index = (n * (n + 1) / 2 + 1);

        schmidtQuasiNorm2 = schmidtQuasiNorm1 * (double)(2 * n - 1) / (double)n;
        schmidtQuasiNorm3 = schmidtQuasiNorm2 * sqrt((double)(n * 2) / (double)(n + 1));
        schmidtQuasiNorm1 = schmidtQuasiNorm2;

        if (n == 1)
        {
            PcupS[n] = PcupS[n - 1];
        }
        else
        {
            k = (double)(((n - 1) * (n - 1)) - 1) / (double)((2 * n - 1) * (2 * n - 3));
            
			PcupS[n] = sin_phi * PcupS[n - 1] - k * PcupS[n - 2];
        }

        /*		  1 nMax  (n+2)    n     m            m           m
                By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                           n=1             m=0   n            n           n  */

        /* Derivative with respect to longitude, divided by radius. */
        By += SphVariables->RelativeRadiusPower[n] *
            (MagneticModel->Secular_Var_Coeff_G[index] * SphVariables->sin_mlambda[1] -
                MagneticModel->Secular_Var_Coeff_H[index] * SphVariables->cos_mlambda[1])
            * PcupS[n] * schmidtQuasiNorm3;
    }

    if (PcupS)
    {
        free(PcupS);
    }
}

void CMagneticResults::RotateMagneticVector(CCoordSpherical* CoordSpherical, CCoordGeodetic* CoordGeodetic, CMagneticResults* MagneticResults)
{
    double Psi;
    
	/* Difference between the spherical and Geodetic latitudes */
    Psi = (M_PI / 180) * (CoordSpherical->phig - CoordGeodetic->phi);

    /* Rotate spherical field components to the Geodetic system */
    Bz = MagneticResults->Bx * sin(Psi) + MagneticResults->Bz * cos(Psi);
    Bx = MagneticResults->Bx * cos(Psi) - MagneticResults->Bz * sin(Psi);
    By = MagneticResults->By;
}

void CMagneticResults::SphericalToCartesian(CCoordSpherical* CoordSpherical, double* x, double* y, double* z)
{
    double radphi;
    double radlambda;

    radphi = CoordSpherical->phig * (M_PI / 180);
    
	radlambda = CoordSpherical->lambda * (M_PI / 180);

    *x = CoordSpherical->r * cos(radphi) * cos(radlambda);
    *y = CoordSpherical->r * cos(radphi) * sin(radlambda);
    *z = CoordSpherical->r * sin(radphi);
}

CGeoMagneticElements::CGeoMagneticElements()
{
    memset(this, 0x00, sizeof(CGeoMagneticElements));
}

CGeoMagneticElements::~CGeoMagneticElements()
{
    if (LegendreFunction)
    {
        delete LegendreFunction;
    }

    if (SphVariables)
    {
        delete SphVariables;
    }

    if (MagneticResultsSph)
    {
        delete MagneticResultsSph;
    }

    if (MagneticResultsSphVar)
    {
        delete MagneticResultsSphVar;
    }

    if (MagneticResultsGeo)
    {
        delete MagneticResultsGeo;
    }

    if (MagneticResultsGeoVar)
    {
        delete MagneticResultsGeoVar;
    }
}

void CGeoMagneticElements::CalculateFieldElements(CEllipsoid* Ellip, CCoordSpherical* CoordSpherical, CCoordGeodetic* CoordGeodetic, CMagneticModel* TimedMagneticModel)
{
    NumTerms = ((TimedMagneticModel->nMax + 1) * (TimedMagneticModel->nMax + 2) / 2);

    LegendreFunction = new CLegendreFunction(NumTerms);

    SphVariables = new CSphericalHarmonicVariables(TimedMagneticModel->nMax);

    SphVariables->ComputeSphericalHarmonicVariables(Ellip, CoordSpherical, TimedMagneticModel->nMax);

    LegendreFunction->AssociatedLegendreFunction(CoordSpherical, TimedMagneticModel->nMax);

    MagneticResultsSph = new CMagneticResults();

    MagneticResultsSph->Summation(LegendreFunction, TimedMagneticModel, SphVariables, CoordSpherical);

    MagneticResultsSphVar = new CMagneticResults();

    MagneticResultsSphVar->SecVarSummation(LegendreFunction, TimedMagneticModel, SphVariables, CoordSpherical);

    MagneticResultsGeo = new CMagneticResults();

    MagneticResultsGeo->RotateMagneticVector(CoordSpherical, CoordGeodetic, MagneticResultsSph);

    MagneticResultsGeoVar = new CMagneticResults();

    MagneticResultsGeoVar->RotateMagneticVector(CoordSpherical, CoordGeodetic, MagneticResultsSphVar);


    CGeoMagneticElements::CalculateGeoMagneticElements(MagneticResultsGeo);

    CGeoMagneticElements::CalculateSecularVariationElements(MagneticResultsGeoVar);
}

void CGeoMagneticElements::CalculateGeoMagneticElements(CMagneticResults* MagneticResultsGeo)
{
    X = MagneticResultsGeo->Bx;
    Y = MagneticResultsGeo->By;
    Z = MagneticResultsGeo->Bz;

    H = sqrt(MagneticResultsGeo->Bx * MagneticResultsGeo->Bx + MagneticResultsGeo->By * MagneticResultsGeo->By);

    F = sqrt(H * H + MagneticResultsGeo->Bz * MagneticResultsGeo->Bz);

    Decl = RAD2DEG(atan2(Y, X));
    Incl = RAD2DEG(atan2(Z, H));
}

void CGeoMagneticElements::CalculateSecularVariationElements(CMagneticResults* MagneticVariation)
{
    Xdot = MagneticVariation->Bx;
    Ydot = MagneticVariation->By;
    Zdot = MagneticVariation->Bz;
    
	Hdot = (X * Xdot + Y * Ydot) / H; 				/* See equation 19 in the WMM technical report */
    Fdot = (X * Xdot + Y * Ydot + Z * Zdot) / F;

    Decldot = 180.0 / M_PI * (X * Ydot - Y * Xdot) / (H * H);
    Incldot = 180.0 / M_PI * (H * Zdot - Z * Hdot) / (F * F);

    GVdot = Decldot;
}

CGeoMagnetic::CGeoMagnetic()
{
    memset(this, 0x00, sizeof(CGeoMagnetic));
}

CGeoMagnetic::~CGeoMagnetic()
{
    delete UserDate;

    delete Geoid;

    delete Ellip;

    delete CoordSpherical;

    delete GeoMagneticElements;

    delete MagneticModel;

    delete TimedMagneticModel;
}

CGeoMagnetic::CGeoMagnetic(double sdate, const char projection)
{
    memset(this, 0x00, sizeof(CGeoMagnetic));

    UserDate = new CDate(sdate);

    Geoid = new CGeoid(projection);

    Ellip = new CEllipsoid();

    CoordSpherical = new CCoordSpherical();

    GeoMagneticElements = new CGeoMagneticElements();

    MagneticModel = new CMagneticModel("WMM.COF");

    TimedMagneticModel = new CMagneticModel(MagneticModel->num_terms);

    TimedMagneticModel->TimelyModifyMagneticModel(UserDate, MagneticModel);
}

void CGeoMagnetic::CalculateFieldElements(double latitude, double longitude, double alt, const char unit)
{
    /* Convert altitude to km */
    switch (unit)
    {
    case 'M':
    {
        alt *= 0.001;

        break;
    }
    
	case 'F':
    {
        alt /= 3280.0839895;

        break;
    }
    }

    CoordGeodetic = new CCoordGeodetic(Geoid, latitude, longitude, alt);

    CoordSpherical->GeodeticToSperical(Ellip, CoordGeodetic);

    GeoMagneticElements->CalculateFieldElements(Ellip, CoordSpherical, CoordGeodetic, TimedMagneticModel);

    delete CoordGeodetic;
}
