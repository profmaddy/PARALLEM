
#include "Data.h"
#include <sys/types.h>
#include <dirent.h>
#include <libgen.h>

int readBurnFDfromFile(Data* data) {

	// Required by gdal n.b. change to c notation
    GDALDatasetH  poDataset;
    GDALRasterBandH  poBand;
    GDALDriverH  poDriver;
    double  adfGeoTransform[6]; //contains projection data
	int  nodatavalue;
	int nXSize, nYSize;
	GDALAllRegister();

    poDataset = GDALOpen( data->FAfile, GA_ReadOnly );
      if( poDataset == NULL )
      {
        printf("File does not exist");
        exit (0);
      }

	poDriver = GDALGetDatasetDriver ( poDataset );

	fprintf(data->outlog, "Driver: %s/%s\n",
			GDALGetDriverShortName( poDriver ),
			GDALGetDriverLongName (poDriver ) );

	poBand = GDALGetRasterBand ( poDataset, 1);
	data->mapInfo.width = GDALGetRasterBandXSize( poBand );
	data->mapInfo.height = GDALGetRasterBandYSize( poBand );

	fprintf(data->outlog, "Size is %d %d\n",
		data->mapInfo.width,
		data->mapInfo.height);

    if( GDALGetGeoTransform( poDataset, adfGeoTransform ) == CE_None )
    {
		data->mapInfo.cellsize = adfGeoTransform[1] ;
		data->mapInfo.xllcorner = adfGeoTransform[0] ;
		data->mapInfo.yllcorner = adfGeoTransform[3] ;
    }

    fprintf(data->outlog,"Cell Size = (%.6f)\n",data->mapInfo.cellsize );

	poBand = GDALGetRasterBand ( poDataset, 1);
	data->mapInfo.nodata = GDALGetRasterNoDataValue(poBand, &nodatavalue);
	fprintf(data->outlog, "No Data Value = (%.6f)\n",data->mapInfo.nodata );

	// Read the values into a buffer zs
	nXSize = GDALGetRasterBandXSize( poBand );
	nYSize = GDALGetRasterBandYSize( poBand );
	// checkCudaErrors(cudaMallocHost((void **)&data->dem, data->mapInfo.height * data->mapInfo.width * sizeof(double)));

	printf("nXSize : %d, nYSize : %d \n", nXSize, nYSize);

    // ************ Read in height data from source ************************
	GDALRasterIO( poBand, GF_Read, 0, 0, nXSize, nYSize, data->fdmod, nXSize, nYSize, GDT_Int32, 0, 0);

    // GDALRasterIO( poBand, GF_Read, 0, 0, nXSize, nYSize, data->fa, nXSize, nYSize, GDT_Float64, 0, 0);
	printf( "SFD Data Read In Complete \n");
	fprintf(data->outlog,"BurnFD Data Read In Complete \n");
	return 1;
}

int readgdalDEMfromFile(Data* data) {

	// Required by gdal n.b. change to c notation
    GDALDatasetH  poDataset;
    GDALRasterBandH  poBand;
    GDALDriverH  poDriver;
    double  adfGeoTransform[6]; //contains projection data
	int  nodatavalue;
	int nXSize, nYSize;

	GDALAllRegister();
    poDataset = GDALOpen( data->demfile, GA_ReadOnly );
      if( poDataset == NULL )
      {
        printf("File does not exist");
        exit (0);
      }

	poDriver = GDALGetDatasetDriver ( poDataset );

	fprintf(data->outlog, "Driver: %s/%s\n",
			GDALGetDriverShortName( poDriver ),
			GDALGetDriverLongName (poDriver ) );

	poBand = GDALGetRasterBand ( poDataset, 1);
	data->mapInfo.width = GDALGetRasterBandXSize( poBand );
	data->mapInfo.height = GDALGetRasterBandYSize( poBand );

	fprintf(data->outlog, "Size is %d %d\n",
		data->mapInfo.width,
		data->mapInfo.height);


    if( GDALGetGeoTransform( poDataset, adfGeoTransform ) == CE_None )
    {
		data->mapInfo.cellsize = adfGeoTransform[1] ;
		data->mapInfo.xllcorner = adfGeoTransform[0] ;
		data->mapInfo.yllcorner = adfGeoTransform[3] ;
    }

    fprintf(data->outlog,"Cell Size = (%.6f)\n",data->mapInfo.cellsize );

	poBand = GDALGetRasterBand ( poDataset, 1);
	data->mapInfo.nodata = GDALGetRasterNoDataValue(poBand, &nodatavalue);
	fprintf(data->outlog, "No Data Value = (%.6f)\n",data->mapInfo.nodata );

	// Read the values into a buffer zs
	nXSize = GDALGetRasterBandXSize( poBand );
	nYSize = GDALGetRasterBandYSize( poBand );
	checkCudaErrors(cudaMallocHost((void **)&data->dem, data->mapInfo.height * data->mapInfo.width * sizeof(double)));

	printf("nXSize : %d, nYSize : %d \n", nXSize, nYSize);

    // ************ Read in height data from source ************************
    GDALRasterIO( poBand, GF_Read, 0, 0, nXSize, nYSize, data->dem, nXSize, nYSize, GDT_Float64, 0, 0);
	printf( "DEM Data Read In Complete \n");
	fprintf(data->outlog,"DEM Data Read In Complete \n");
	return 1;
}

int readgdalBedrockfromFile(Data* data, Data* device) {

	// Required by gdal n.b. change to c notation
    GDALDatasetH  poDataset;
    GDALRasterBandH  poBand;
    GDALDriverH  poDriver;
    int nXSize, nYSize;

	GDALAllRegister();
    poDataset = GDALOpen( data->bedrockfile, GA_ReadOnly );
      if( poDataset == NULL )
      {
        printf("Bedrock file does not exist");
        exit (0);
      }
	poDriver = GDALGetDatasetDriver ( poDataset );
   	poBand = GDALGetRasterBand ( poDataset, 1);
	nXSize = GDALGetRasterBandXSize( poBand );
	nYSize = GDALGetRasterBandYSize( poBand );

	//checkCudaErrors(cudaMallocHost((void **)&data->Flow_C, nXSize * nYSize * sizeof(double)));

	// changed to load thickness data on 3/10/18

	checkCudaErrors(cudaMallocHost((void **)&data->soilTPtr, nXSize * nYSize * sizeof(double)));

	    // ************ Read in bedrock data from source ************************
    GDALRasterIO( poBand, GF_Read, 0, 0, nXSize, nYSize, data->soilTPtr, nXSize, nYSize, GDT_Float64, 0, 0 );
	printf( "Soil Thickness Data Read In Complete \n");
	fprintf(data->outlog,"Soil Thickness Data Read In Complete \n");

	//checkCudaErrors( cudaMalloc( (void**) &device->soilTPtr, nXSize * nYSize * sizeof(double)) );
	//checkCudaErrors( cudaMemcpy ( device->soilTPtr,  data->soilTPtr,  nXSize * nYSize * sizeof(double), cudaMemcpyHostToDevice) );
	//printf( "Soil Thickness Data sent to GPU \n");
	//fprintf(data->outlog,"Bedrock Data sent to GPU \n");

	return 1;
}

int readgdalmaskfromFile(Data* data) {

	// Required by gdal n.b. change to c notation
    GDALDatasetH  poDataset;
    GDALRasterBandH  poBand;
    GDALDriverH  poDriver;
    int  nodatavalue;
	int nXSize, nYSize;
	GDALAllRegister();

    poDataset = GDALOpen( data->maskfile, GA_ReadOnly );
      if( poDataset == NULL )
      {
        printf("File does not exist");
        exit(0);
      }

	poDriver = GDALGetDatasetDriver ( poDataset );

	fprintf(data->outlog,"Driver: %s/%s\n",
			GDALGetDriverShortName( poDriver ),
			GDALGetDriverLongName (poDriver ) );

	poBand = GDALGetRasterBand ( poDataset, 1);
	nodatavalue = GDALGetRasterNoDataValue(poBand, &nodatavalue);
	printf( "No Data Value = (%d)\n",nodatavalue );

	// Read the values into a buffer zs
	nXSize = GDALGetRasterBandXSize( poBand );
	nYSize = GDALGetRasterBandYSize( poBand );
	data->mask = (int *) CPLMalloc(sizeof(int)*nXSize*nYSize);

	fprintf(data->outlog,"nXSize : %d, nYSize : %d \n", nXSize, nYSize);

    // ************ Read in height data from source ************************
    GDALRasterIO( poBand, GF_Read, 0, 0, nXSize, nYSize, data->mask, nXSize, nYSize, GDT_Int32, 0, 0 );
	printf( "Mask Data Read In Complete \n");
	fprintf(data->outlog,"Mask Data Read In Complete \n");
	return 1;
}

void savegrids(Data* data, Catchment* catchments,  int iteration)
{
	int i = iteration;
	if (data->flowdirectiontype == 0) sprintf(data->FDfile, "%s/%d_SFD.tif",data->matrixDIR,  i);
	if (data->flowdirectiontype == 1) sprintf(data->FDfile, "%s/%d_MFD.tif",  data->matrixDIR, i);
	writegdalGRIDtoFile(data, catchments, data->FDfile, 1, 0);

	// Write the recalculated FA
	sprintf(data->FAfile, "%s/%d_FA.tif", data->matrixDIR,  i);
	writegdalGRIDtoFile(data, catchments, data->FAfile, 0, 1);


	if ( (i>=2) & (i%20==0) ) {

	budget_calc(data);
	// Write the recalculated DEM
	sprintf(data->heightfile, "%s/%d_height.tif", data->matrixDIR, i);
	writegdalGRIDtoFile(data, catchments, data->heightfile, 0, 0);

	// Write the budget matrix
	sprintf(data->diff_file, "%s/%d_diff.tif", data->matrixDIR,  i);
	writegdalGRIDtoFile(data, catchments, data->diff_file, 0, 18);

	copyprevht(data);

	}

/**************************************************************/
	if ( (i>=2) & (i%500==0) ) {

		// Write the recalculated SFD
		if (data->flowdirectiontype == 0) sprintf(data->FDfile, "%s/%d_SFD.tif",data->matrixDIR,  i);
		if (data->flowdirectiontype == 1) sprintf(data->FDfile, "%s/%d_MFD.tif",  data->matrixDIR, i);
		writegdalGRIDtoFile(data, catchments, data->FDfile, 1, 0);

		// Write the recalculated FA
		sprintf(data->FAfile, "%s/%d_FA.tif", data->matrixDIR,  i);
		writegdalGRIDtoFile(data, catchments, data->FAfile, 0, 1);

		// Write the recalculated erosion
		sprintf(data->erofile, "%s/%d_ero.tif", data->matrixDIR,   i);
		writegdalGRIDtoFile(data, catchments, data->erofile, 0, 2);

		// Write the recalculated incision
		sprintf(data->incifile, "%s/%d_inci.tif", data->matrixDIR, i);
		writegdalGRIDtoFile(data, catchments, data->incifile, 0, 3);

		// Write the recalculated gelifluction
		sprintf(data->gelifile, "%s/%d_geli.tif", data->matrixDIR,   i);
		writegdalGRIDtoFile(data, catchments, data->gelifile, 0, 13);

		// Write the recalculated deposition
		sprintf(data->depofile, "%s/%d_depo.tif", data->matrixDIR, i);
		writegdalGRIDtoFile(data, catchments, data->depofile, 0, 4);

		// Write the recalculated slope
		sprintf(data->slopefile, "%s/%d_slope.tif", data->matrixDIR,   i);
		writegdalGRIDtoFile(data, catchments, data->slopefile, 0, 5);

		sprintf(data->Precipfile, "%s/%d_precip.tif", data->matrixDIR, iteration);
		writegdalGRIDtoFile(data, catchments, data->Precipfile, 0, 16);

		sprintf(data->Tempfile, "%s/%d_temp.tif", data->matrixDIR, iteration);
		writegdalGRIDtoFile(data, catchments, data->Tempfile, 0, 17);

		// Write the recalculated soil Thickness values
		sprintf(data->soilTfile, "%s/%d_soilT.tif", data->matrixDIR, i);
		writegdalGRIDtoFile(data, catchments, data->soilTfile, 0, 9);

		// Write the recalculated stone
		sprintf(data->stonesfile, "%s/%d_stone.tif", data->matrixDIR,  i);
		writegdalGRIDtoFile(data, catchments, data->stonesfile, 0, 7);

		// Write the recalculated fines
		sprintf(data->finesfile, "%s/%d_fines.tif", data->matrixDIR,   i);
		writegdalGRIDtoFile(data, catchments, data->finesfile, 0, 6);

	}
	/**************************************************************/

	if ( (i>=2) & (i%1000==0) )
	{
			// Write the recalculated soil Moisture values
			sprintf(data->soilMfile, "%s/%d_soilM.tif", data->matrixDIR,   i);
			writegdalGRIDtoFile(data, catchments, data->soilMfile, 0, 15);

			// Write the recalculated flatmask
			sprintf(data->flatfile, "%s/%d_flatmask.tif", data->matrixDIR, i);
			writegdalGRIDtoFile(data, catchments, data->flatfile, 1, 4);

			// Write the recalculated soil Bio values
			sprintf(data->soilBfile, "%s/%d_soilB.tif", data->matrixDIR,   i);
			writegdalGRIDtoFile(data, catchments, data->soilBfile, 0, 14);

			// Write the recalculated soil nutrients values
			sprintf(data->nutfile, "%s/%d_nut.tif",data->matrixDIR,  i);
			writegdalGRIDtoFile(data, catchments, data->nutfile, 0, 10);

			// Write the recalculated total Bio  values
			sprintf(data->totbiofile, "%s/%d_totbio.tif", data->matrixDIR,   i);
			writegdalGRIDtoFile(data, catchments, data->totbiofile, 0, 8);

			// Write the recalculated chemical weathering
			sprintf(data->wCfile, "%s/%d_wC.tif", data->matrixDIR,  i);
			writegdalGRIDtoFile(data, catchments, data->wCfile, 0, 11);

			// Write the recalculated physical weathering
			sprintf(data->wPfile, "%s/%d_wP.tif",  data->matrixDIR, i);
			writegdalGRIDtoFile(data, catchments, data->wPfile, 0, 12);

			//sprintf(contrib, "output1/contrib_%dm_%d.tif", data.WhichGrid, i);
			//writegdalGRIDtoFile(&data, &catchments, contrib, 1, 3);
    }

	return;
}

//int writegdalGRIDtoFile(Data* data, const char* fileName, int whichtype, int what) {
int writeGRIDtoFile(Data* data, char* fileName, int whichtype, int what) {

	double transformcell = - (data->mapInfo.cellsize) ;
	double adfGeoTransform[6] = {data->mapInfo.xllcorner, data->mapInfo.cellsize, 0, data->mapInfo.yllcorner, 0, transformcell };


	//printf(" %f, %f, %f, %f, %f, %f ", adfGeoTransform[0],adfGeoTransform[1],adfGeoTransform[2],adfGeoTransform[3],adfGeoTransform[4],adfGeoTransform[5]);

	GDALDatasetH DDstDS; // for doubles
	GDALDatasetH IntDstDS; // for integers

    OGRSpatialReferenceH hSRS;
    char *pszSRS_WKT = NULL;

    char **Options = NULL;
	const char *pszFormat = "GTiff";
    GDALDriverH hDriver = GDALGetDriverByName( pszFormat );
	GDALRasterBandH  DstBand, DstBandFD;

	if( hDriver == NULL )  exit( 1 );

	Options = CSLSetNameValue( Options, "COMPRESS","LZW");

	int espgcode = 25830 ; // this is ETRS89 / UTM Zon 30N from 27700 for OS GB grid;
    hSRS = OSRNewSpatialReference( NULL );
    OSRSetWellKnownGeogCS( hSRS, "EPSG:27700" );
    OSRImportFromEPSG(hSRS, espgcode);
    OSRExportToWkt( hSRS, &pszSRS_WKT );
    OSRDestroySpatialReference( hSRS );

	switch (whichtype) {
		case 0: // write doubles

		DDstDS = GDALCreate( hDriver, fileName, data->mapInfo.width, data->mapInfo.height, 1, GDT_Float64, Options );
		GDALSetGeoTransform( DDstDS, adfGeoTransform );


	    GDALSetProjection( DDstDS, pszSRS_WKT );
	    //CPLFree( pszSRS_WKT );

		DstBand = GDALGetRasterBand( DDstDS, 1 );
		GDALSetRasterNoDataValue(DstBand, -9999);


			switch (what) {
				case 0: // write the DEM
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->dem, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 1: // write FA
					GDALSetRasterNoDataValue(DstBand, -9999);
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->fa, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 2: // write erosion (diffuse)
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->eroPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 3: // write incision (con concentrated)
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->inciPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 4: // write deposition
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->depoPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 5: // write slope
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->SlopePtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 6: // write fines
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->finesPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 7: // write stone
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->stonePtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 8: // write total Bio
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->TotBPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 9: // write soil thickness
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->soilTPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 10: // write nutrients
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->nutPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 11: // write chemical weathering
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->weatherC, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 12: // write Physical weathering
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->weatherP, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 13: // write solifluction
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->geliPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 14: // write soilBio
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->soilBPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 15: // write soil Moisture
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->soilMPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 16: // write Precipfile
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->rainmat, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 17: // write Tempfile
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->tempmat, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 18: // write dzfile
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->shortest_paths, data->mapInfo.width, data->mapInfo.height, GDT_Float32, 0, 0 );
				break;
			}


		GDALClose (DDstDS);

		break;

		case 1: // write the integers
		IntDstDS = GDALCreate( hDriver, fileName, data->mapInfo.width, data->mapInfo.height, 1, GDT_Int32, Options );
		printf("file name for FD file is %s\n", fileName);
		GDALSetGeoTransform( IntDstDS, adfGeoTransform );
	    GDALSetProjection( IntDstDS, pszSRS_WKT );
	    CPLFree( pszSRS_WKT );


		DstBandFD = GDALGetRasterBand( IntDstDS, 1 );
			switch (what) {
			case 0: // write the FD file
				    GDALSetRasterNoDataValue(DstBandFD, -9999);
					GDALRasterIO( DstBandFD, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->SFD, data->mapInfo.width, data->mapInfo.height, GDT_Int32, 0, 0 );
			break;

			case 1: // write the watershed file
				GDALSetRasterNoDataValue(DstBandFD, -9999);
				GDALRasterIO( DstBandFD, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->watershed_id, data->mapInfo.width, data->mapInfo.height, GDT_Int32, 0, 0 );
				// write the catchment id map
				//GDALRasterIO( DstBandFD, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, catchments->watershed_id, data->mapInfo.width, data->mapInfo.height, GDT_Int32, 0, 0 );
			break;

			case 2: // write the catchment mask
				//GDALRasterIO( DstBandFD, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, catchments->mask, data->mapInfo.width, data->mapInfo.height, GDT_Int32, 0, 0 );
			break;

			case 3: // write the FA okgrid
				GDALRasterIO( DstBandFD, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->contribA, data->mapInfo.width, data->mapInfo.height, GDT_Int32, 0, 0 );
			break;

			case 4: // write the flatmask
				GDALRasterIO( DstBandFD, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->flatmask, data->mapInfo.width, data->mapInfo.height, GDT_Int32, 0, 0 );
			break;
			}
		GDALClose (IntDstDS);
		break;
	}

	return 1;
}

//int writegdalGRIDtoFile(Data* data, const char* fileName, int whichtype, int what) {
int writegdalGRIDtoFile(Data* data, Catchment* catchments, char* fileName, int whichtype, int what) {

	double transformcell = - (data->mapInfo.cellsize) ;
	double adfGeoTransform[6] = {data->mapInfo.xllcorner, data->mapInfo.cellsize, 0, data->mapInfo.yllcorner, 0, transformcell };


	//printf(" %f, %f, %f, %f, %f, %f ", adfGeoTransform[0],adfGeoTransform[1],adfGeoTransform[2],adfGeoTransform[3],adfGeoTransform[4],adfGeoTransform[5]);

	GDALDatasetH DDstDS; // for doubles
	GDALDatasetH IntDstDS; // for integers

    OGRSpatialReferenceH hSRS;
    char *pszSRS_WKT = NULL;

    char **Options = NULL;
	const char *pszFormat = "GTiff";
    GDALDriverH hDriver = GDALGetDriverByName( pszFormat );
	GDALRasterBandH  DstBand, DstBandFD;

	if( hDriver == NULL )  exit( 1 );

	Options = CSLSetNameValue( Options, "COMPRESS","LZW");

	int espgcode = 25830 ; // this is ETRS89 / UTM Zon 30N from 27700 for OS GB grid;
    hSRS = OSRNewSpatialReference( NULL );
    //OSRSetWellKnownGeogCS( hSRS, "EPSG:27700" );
    OSRImportFromEPSG(hSRS, espgcode);
    OSRExportToWkt( hSRS, &pszSRS_WKT );
    OSRDestroySpatialReference( hSRS );

	switch (whichtype) {
		case 0: // write doubles
		DDstDS = GDALCreate( hDriver, fileName, data->mapInfo.width, data->mapInfo.height, 1, GDT_Float64, Options );
		GDALSetGeoTransform( DDstDS, adfGeoTransform );


	    GDALSetProjection( DDstDS, pszSRS_WKT );
	    //CPLFree( pszSRS_WKT );

		DstBand = GDALGetRasterBand( DDstDS, 1 );
		GDALSetRasterNoDataValue(DstBand, -9999);


			switch (what) {
				case 0: // write the DEM
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->dem, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 1: // write FA
					GDALSetRasterNoDataValue(DstBand, -9999);
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->fa, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 2: // write erosion (diffuse)
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->eroPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 3: // write incision (con concentrated)
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->inciPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 4: // write deposition
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->depoPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 5: // write slope
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->SlopePtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 6: // write fines
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->finesPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 7: // write stone
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->stonePtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 8: // write total Bio
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->TotBPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 9: // write soil thickness
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->soilTPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 10: // write nutrients
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->nutPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 11: // write chemical weathering
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->weatherC, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 12: // write Physical weathering
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->weatherP, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 13: // write solifluction
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->geliPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 14: // write soilBio
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->soilBPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 15: // write soil Moisture
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->soilMPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 16: // write Precipfile
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->rainmat, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 17: // write Tempfile
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->tempmat, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 18: // write dzfile
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->diffdem, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;
			}


		GDALClose (DDstDS);

		break;

		case 1: // write the integers
		IntDstDS = GDALCreate( hDriver, fileName, data->mapInfo.width, data->mapInfo.height, 1, GDT_Int32, Options );
		GDALSetGeoTransform( IntDstDS, adfGeoTransform );

	    GDALSetProjection( IntDstDS, pszSRS_WKT );
	    CPLFree( pszSRS_WKT );


		DstBandFD = GDALGetRasterBand( IntDstDS, 1 );
			switch (what) {
			case 0: // write the FD file
				    GDALSetRasterNoDataValue(DstBandFD, -9999);
					GDALRasterIO( DstBandFD, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->fd, data->mapInfo.width, data->mapInfo.height, GDT_Int32, 0, 0 );
			break;

			case 1: // write the catchment id map
				GDALRasterIO( DstBandFD, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, catchments->watershed_id, data->mapInfo.width, data->mapInfo.height, GDT_Int32, 0, 0 );
			break;

			case 2: // write the catchment mask
				GDALRasterIO( DstBandFD, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, catchments->mask, data->mapInfo.width, data->mapInfo.height, GDT_Int32, 0, 0 );
			break;

			case 3: // write the contributing area
				GDALRasterIO( DstBandFD, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->contribA, data->mapInfo.width, data->mapInfo.height, GDT_Int32, 0, 0 );
			break;

			case 4: // write the flatmask
				GDALRasterIO( DstBandFD, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->flatmask, data->mapInfo.width, data->mapInfo.height, GDT_Int32, 0, 0 );
			break;
			}
		GDALClose (IntDstDS);
		break;
	}

	return 1;
}


int retriveProcessMatrices(Data* data) //override matrices with existing data

{
	// Required by gdal n.b. change to c notation
    GDALDatasetH  height, stone, fines, soilM, soilB, soilT, Nut, TotB ;
    GDALRasterBandH  heightB, stoneB, finesB, soilMB, soilBB, soilTB, NutB, TotBB;
    GDALDriverH  poDriver;
    double  adfGeoTransform[6]; //contains projection data
	int  nodatavalue;
	int nXSize, nYSize;


	char* fileName = (char*) malloc(sizeof(char) *100);

	nXSize = data->mapInfo.width;
	nYSize = data->mapInfo.height;

	GDALAllRegister();

	// ************ Read in height data  from file ************************

	    sprintf(fileName, "restart/%d_height.tif", data->start_iter);
	    height = GDALOpen( fileName, GA_ReadOnly );
	      if( height == NULL )
	      {
	        printf("height File does not exist");
	        exit (0);
	      }
		poDriver = GDALGetDatasetDriver ( height ); // this only needs to be done once as all files in same format

		fprintf(data->outlog, "Driver: %s/%s\n",
				GDALGetDriverShortName( poDriver ),
				GDALGetDriverLongName (poDriver ) );

	heightB = GDALGetRasterBand ( height, 1);
    GDALRasterIO( heightB, GF_Read, 0, 0, nXSize, nYSize, data->dem, nXSize, nYSize, GDT_Float64  , 0, 0 );
	printf( "Height Data Read In Complete \n");
	fprintf(data->outlog,"Height Data Read In Complete \n");


    // ************ Read in stone data  from file ************************

	sprintf(fileName, "restart/%d_stone.tif", data->start_iter);
	stone = GDALOpen( fileName, GA_ReadOnly );
	      if( stone == NULL )
	      {
	        printf("stone File does not exist");
	        exit (0);
	      }
	stoneB = GDALGetRasterBand ( stone, 1);
    GDALRasterIO( stoneB, GF_Read, 0, 0, nXSize, nYSize, data->stonePtr, nXSize, nYSize, GDT_Float64  , 0, 0 );
	printf( "Stone Data Read In Complete \n");
	fprintf(data->outlog,"Stone Data Read In Complete \n");

    // ************ Read in fines data  from file ************************

    sprintf(fileName, "restart/%d_fines.tif", data->start_iter);
    fines = GDALOpen( fileName, GA_ReadOnly );
      if( fines == NULL )
      {
        printf("fines File does not exist");
        exit (0);
      }
	finesB = GDALGetRasterBand ( fines, 1);
	GDALRasterIO( finesB, GF_Read, 0, 0, nXSize, nYSize, data->finesPtr, nXSize, nYSize, GDT_Float64  , 0, 0 );
	printf( "Fines Data Read In Complete \n");
	fprintf(data->outlog,"Fines Data Read In Complete \n");

    // ************ Read in Soil Moisturedata  from file ************************

    sprintf(fileName, "restart/%d_soilM.tif", data->start_iter);
    soilM = GDALOpen( fileName, GA_ReadOnly );
      if( soilM == NULL )
      {
        printf("soilM File does not exist");
        exit (0);
      }
	soilMB = GDALGetRasterBand ( soilM, 1);
	GDALRasterIO( soilMB, GF_Read, 0, 0, nXSize, nYSize, data->soilMPtr, nXSize, nYSize, GDT_Float64  , 0, 0 );
	printf( "Soil Moisture Data Read In Complete \n");
	fprintf(data->outlog,"Soil Moisture Data Read In Complete \n");

    // ************ Read in Soil Biodata  from file ************************

    sprintf(fileName, "restart/%d_soilB.tif", data->start_iter);
    soilB = GDALOpen( fileName, GA_ReadOnly );
      if( soilB == NULL )
      {
        printf("soilB File does not exist");
        exit (0);
      }
	soilBB = GDALGetRasterBand ( soilB, 1);
	GDALRasterIO( soilBB, GF_Read, 0, 0, nXSize, nYSize, data->soilBPtr, nXSize, nYSize, GDT_Float64  , 0, 0 );
	printf( "Soil Bio Data Read In Complete \n");
	fprintf(data->outlog,"Soil Bio Data Read In Complete \n");

    // ************ Read in Soil Thickness data  from file ************************

    sprintf(fileName, "restart/%d_soilT.tif", data->start_iter);
    soilT = GDALOpen( fileName, GA_ReadOnly );
      if( soilT == NULL )
      {
        printf("soilT File does not exist");
        exit (0);
      }
	soilTB = GDALGetRasterBand ( soilT, 1);
	GDALRasterIO( soilTB, GF_Read, 0, 0, nXSize, nYSize, data->soilTPtr, nXSize, nYSize, GDT_Float64  , 0, 0 );
	printf( "Soil Thickness Data Read In Complete \n");
	fprintf(data->outlog,"Soil Thickness Data Read In Complete \n");

    // ************ Read in Nutrient data  from file ************************
    sprintf(fileName, "restart/%d_nut.tif", data->start_iter);
    Nut = GDALOpen( fileName, GA_ReadOnly );
      if( Nut == NULL )
      {
        printf("Nut File does not exist");
        exit (0);
      }
	NutB = GDALGetRasterBand ( Nut, 1);
	GDALRasterIO( NutB, GF_Read, 0, 0, nXSize, nYSize, data->nutPtr, nXSize, nYSize, GDT_Float64  , 0, 0 );
	printf( "Nutrient Data Read In Complete \n");
	fprintf(data->outlog,"Nutrient Data Read In Complete \n");

    // ************ Read in Total Bio data  from file ************************

    sprintf(fileName, "restart/%d_totbio.tif", data->start_iter);
    TotB = GDALOpen( fileName, GA_ReadOnly );
      if( TotB == NULL )
      {
        printf("TotB File does not exist");
        exit (0);
      }
	TotBB = GDALGetRasterBand (TotB, 1);
	GDALRasterIO( TotBB, GF_Read, 0, 0, nXSize, nYSize, data->TotBPtr, nXSize, nYSize, GDT_Float64  , 0, 0 );
	printf( "Total Bio Data Read In Complete \n");
	fprintf(data->outlog,"Total Bio Data Read In Complete \n");

	printf( "\n\nAll restart Data Read In\n\n");

	data->runiniter = 0;

	return 1;
}
