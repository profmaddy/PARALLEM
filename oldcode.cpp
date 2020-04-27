void DirFlow(Data* data)
{
    int cellx, celly;
    int dcell;
    int aspect, aspectsfd, newaspect;

    double thiscellht, targetcellht;
    int upslopecount;
    int downslopecount;
    int flatcount;
    int flatcounter;
    int sinkcounter;
    int mycounts;

    double smax, stemp, slopetot;
    float dc;

    stemp = 0;
    smax = 0.0;// needed for SFD aspect
    flatcounter = 0;
    sinkcounter = 0;

    double cell_size;
    cell_size = data->mapInfo.cellsize;
    int ncell_y = data->mapInfo.height;
    int ncell_x = data->mapInfo.width;
    int self;

    int xmove[9] = { 0,1,1,1,0,-1,-1,-1,0 };
    int ymove[9] = { -1,-1,0,1,1,1,0,-1,0 };
    int* dx;
    int* dy;
    dx = &xmove[0];
    dy = &ymove[0];

    data->dx = dx;
    data->dy = dy;

    for (int irow = 0; irow < ncell_y; irow++)  //irow loop
    {
        for (int icol = 0; icol < ncell_x; icol++) //icol loop
        {
            flatcount = 0;
            slopetot = 0.0; // for down slope directions and mfd only
            aspect = 0; // we will now code using MFD aspect coding -default aspect is none i.e. zero, this will signify sink or flat for routing
            aspectsfd = 0;
            upslopecount = 0;
            downslopecount = 0;
            self = irow * ncell_x + icol;
            //printf("self %d\n", self);

            if (data->dem[self] > 0) // in dem loop
            {
                for (dcell = 0; dcell < 8; dcell++)
                { //dcell loop
                    cellx = icol + data->dx[dcell];
                    celly = irow + data->dy[dcell];
                    if (cellx >= 0 && cellx < ncell_x && celly >= 0 && celly < ncell_y) // for each of my neighbours
                    {
                        newaspect = 2 ^ ((dcell + 6) % 8);
                        /*if (dcell == 0)  newaspect = 64;
                         if (dcell == 1)  newaspect = 128;
                         if (dcell == 2)  newaspect = 1;
                         if (dcell == 3)  newaspect = 2;
                         if (dcell == 4)  newaspect = 4;
                         if (dcell == 5)  newaspect = 8;
                         if (dcell == 6)  newaspect = 16;
                         if (dcell == 7)  newaspect = 32;*/
                        if (dcell % 2 == 0)
                            dc = 1.0;   //even directions are cardinal
                        else
                            dc = 1.41;  //odd directions are diagonals
                      //calculate the slope
                        thiscellht = data->dem[self];
                        targetcellht = data->dem[celly * ncell_x + cellx];
                        printf("thiscell %f, targetcell %f \n", thiscellht, targetcellht);
                        stemp = (thiscellht - targetcellht) / (cell_size * dc);
                        //printf("calculated slope = %f \n", stemp);
                        if (targetcellht == -9999)
                        {
                            stemp = -0.0000001;
                        } // needed to identify sinks on boundary

//store the slope and aspect value if it flows here i.e. down hill
                        if (thiscellht > targetcellht)
                        {
                            downslopecount += 1;
                            data->Slopes[self + dcell] = stemp;
                            slopetot += data->Slopes[self + dcell];
                            if (stemp > smax) // default aspect set to aspect of max slope for SFD
                            {
                                smax = stemp;
                                aspectsfd = newaspect;
                                //if (aspectsfd == 0) return;
                            }
                            aspect = aspect + newaspect; // for MFD

                            //if (aspect == 0) return;
                        }
                        if (thiscellht < targetcellht) upslopecount += 1;
                        if (thiscellht == targetcellht)   flatcount += 1;
                    }
                }
            } //dcell loop i.e. look at my neighbours

            if (upslopecount > 7) {
                sinkcounter = sinkcounter + 1;
                aspect = 0;
                aspectsfd = 0;
            } // I got nowhere to go because I am a sink

            if ((flatcount + upslopecount) > 7)
            {
                if (upslopecount <= 7) { // do not double count the sinks
                    flatcounter = flatcounter + 1;
                    aspect = 0;
                    aspectsfd = 0;
                }
                // I got nowhere to go because I am a flat next to either adjacent higher or equal height cells i.e. no lower cells

                if (downslopecount > 0) {
                    printf("upslope %d, downslope %d, flatcounter %d aspect %d \n", upslopecount, downslopecount, flatcount, aspect);
                }


                mycounts = flatcount + upslopecount + downslopecount;
                if (mycounts != 8) printf("something went wrong, mycounts = %d", mycounts);

                data->fd[self] = aspect;
                data->SFD[self] = aspectsfd;// set the SFD direction - this is used for all processes i.e. they are not MFD

                flatcount = 0;
                slopetot = 0.0; // for down slope directions and mfd only
                aspect = 0; // we will now code using MFD aspect coding -default aspect is none i.e. zero, this will signify sink or flat for routing
                aspectsfd = 0;
                upslopecount = 0;
                downslopecount = 0;

            } // in dem loop
        } // icol
    }// irow

    printf("Sinkcounter = %d \n", sinkcounter);
    printf("flatcounter = %d \n", flatcounter);
}

void checkslopeandprop(Data* data)
{
    int self;
    int idx;
    int counter;

    counter = 0;

    for (int irow = 0; irow < data->mapInfo.height; irow++)  //irow loop
    {
        for (int icol = 0; icol < data->mapInfo.width; icol++) //icol loop
        {
            self = self = irow * data->mapInfo.height + icol;

            if (data->mask[self] == 1)
            {

                for (int dcell = 0; dcell < 8; dcell++)
                { //dcell loop
                    if ((data->prop[self + dcell] > 0) && (data->Slopes[self + dcell] == 0.0))
                        printf(" prop allocated to direction with no slope %d \n", self);

                    idx = self + dcell;

                    if ((data->fd[self] == 1) && (data->prop[idx] != 1)) {
                        printf("SFD %d has proportion %lf \n", data->fd[self], data->prop[idx]);
                        counter++;
                    }

                    if ((data->fd[self] == 2) && (data->prop[idx] != 1)) {
                        printf("SFD %d has proportion %lf \n", data->fd[self], data->prop[idx]);
                        counter++;
                    }

                    if ((data->fd[self] == 4) && (data->prop[idx] != 1)) {
                        printf("SFD %d has proportion %lf \n", data->fd[self], data->prop[idx]);
                        counter++;
                    }

                    if ((data->fd[self] == 8) && (data->prop[idx] != 1)) {
                        printf("SFD %d has proportion %lf \n", data->fd[self], data->prop[idx]);
                        counter++;
                    }

                    if ((data->fd[self] == 16) && (data->prop[idx] != 1)) {
                        printf("SFD %d has proportion %lf \n", data->fd[self], data->prop[idx]);
                        counter++;
                    }

                    if ((data->fd[self] == 32) && (data->prop[idx] != 1)) {
                        printf("SFD %d has proportion %lf \n", data->fd[self], data->prop[idx]);
                        counter++;
                    }

                    if ((data->fd[self] == 64) && (data->prop[idx] != 1)) {
                        printf("SFD %d has proportion %lf \n", data->fd[self], data->prop[idx]);
                        counter++;
                    }

                    if ((data->fd[self] == 128) && (data->prop[idx] != 1)) {
                        printf("SFD %d has proportion %lf \n", data->fd[self], data->prop[idx]);
                        counter++;
                    }

                }
            }
        }
    }
    printf("prop and slope check complete %d \n", counter);
}

__device__ double mfdsed(int client, int self, double* hv, int selffd, int gridCols)
{
    int do001, do002, do004, do008, do016, do032, do064, do128;
    double slope;
    double slopetotal;
    double slope001, slope002, slope004, slope008, slope016, slope032, slope064, slope128;
    double prop; // proportion
    int targets;
    double selfhv; //hv = height value
    double thedistance;

    thedistance = 1.0;
    targets = 0;
    slopetotal = 0.0;
    selfhv = hv[self];

    // filter ONLY the directions indicated in the MFD i.e. downslope directions and set 1 otherwise 0
    do001 = (selffd & 1) != 0;
    do002 = (selffd & 2) != 0;
    do004 = (selffd & 4) != 0;
    do008 = (selffd & 8) != 0;
    do016 = (selffd & 16) != 0;
    do032 = (selffd & 32) != 0;
    do064 = (selffd & 64) != 0;
    do128 = (selffd & 128) != 0;

    // all slope values need to be positive for sum to proportion to work
    slope001 = do001 * abs((double)selfhv - hv[self + 1]); // calculate the slopes
    slope002 = do002 * abs((double)selfhv - hv[self + gridCols + 1]) / 1.41; // note the change in denominator for diagonals.
    slope004 = do004 * abs((double)selfhv - hv[self + gridCols]);
    slope008 = do008 * abs((double)selfhv - hv[self + gridCols - 1]) / 1.41;
    slope016 = do016 * abs((double)selfhv - hv[self - 1]);
    slope032 = do032 * abs((double)selfhv - hv[self - gridCols - 1]) / 1.41;
    slope064 = do064 * abs((double)selfhv - hv[self - gridCols]);
    slope128 = do128 * abs((double)selfhv - hv[self - gridCols + 1]) / 1.41;

    if ((do001 == 1) && (slope001 == 0)) slope001 = 0.01;
    if ((do002 == 1) && (slope002 == 0)) slope002 = 0.01;
    if ((do004 == 1) && (slope004 == 0)) slope004 = 0.01;
    if ((do008 == 1) && (slope008 == 0)) slope008 = 0.01;
    if ((do016 == 1) && (slope016 == 0)) slope016 = 0.01;
    if ((do032 == 1) && (slope032 == 0)) slope032 = 0.01;
    if ((do064 == 1) && (slope064 == 0)) slope064 = 0.01;
    if ((do128 == 1) && (slope128 == 0)) slope128 = 0.01;

    slopetotal = (slope001 + slope004 + slope016 + slope064) + (slope002 + slope008 + slope032 + slope128);  //diagonals have already been divided by 1.4 see above

    targets = do001 + do002 + do004 + do008 + do016 + do032 + do064 + do128;

    if (do002 == 1 || do008 == 1 || do032 == 1 || do128 == 1) thedistance = 1.41;

    // when pulling slope would be negative - need to make sure it is positive for proportions
    slope = ((double)selfhv - hv[client]) / thedistance; // placed back in Aug 30th 17

    //if ((slope<=0)|| (slope>6)) slope = 6;
    if ((slope == 0)) slope = 0.01;

    prop = 1.0 / targets; // default for SFD

 /* if ( targets > 1)
    {
    prop = (double) (slope/(slopetotal));
    } */

    return prop;
}
__device__ double mfdsed(int client, int self, double* hv, int selffd, int gridCols)
{
    int do001, do002, do004, do008, do016, do032, do064, do128;
    double slope;
    double slopetotal;
    double slope001, slope002, slope004, slope008, slope016, slope032, slope064, slope128;
    double prop; // proportion
    int targets;
    double selfhv; //hv = height value
    double thedistance;

    thedistance = 1.0;
    targets = 0;
    slopetotal = 0.0;
    selfhv = hv[self];

    // filter ONLY the directions indicated in the MFD i.e. downslope directions and set 1 otherwise 0
    do001 = (selffd & 1) != 0;
    do002 = (selffd & 2) != 0;
    do004 = (selffd & 4) != 0;
    do008 = (selffd & 8) != 0;
    do016 = (selffd & 16) != 0;
    do032 = (selffd & 32) != 0;
    do064 = (selffd & 64) != 0;
    do128 = (selffd & 128) != 0;

    // all slope values need to be positive for sum to proportion to work
    slope001 = do001 * abs((double)selfhv - hv[self + 1]); // calculate the slopes
    slope002 = do002 * abs((double)selfhv - hv[self + gridCols + 1]) / 1.41; // note the change in denominator for diagonals.
    slope004 = do004 * abs((double)selfhv - hv[self + gridCols]);
    slope008 = do008 * abs((double)selfhv - hv[self + gridCols - 1]) / 1.41;
    slope016 = do016 * abs((double)selfhv - hv[self - 1]);
    slope032 = do032 * abs((double)selfhv - hv[self - gridCols - 1]) / 1.41;
    slope064 = do064 * abs((double)selfhv - hv[self - gridCols]);
    slope128 = do128 * abs((double)selfhv - hv[self - gridCols + 1]) / 1.41;

    if ((do001 == 1) && (slope001 == 0)) slope001 = 0.01;
    if ((do002 == 1) && (slope002 == 0)) slope002 = 0.01;
    if ((do004 == 1) && (slope004 == 0)) slope004 = 0.01;
    if ((do008 == 1) && (slope008 == 0)) slope008 = 0.01;
    if ((do016 == 1) && (slope016 == 0)) slope016 = 0.01;
    if ((do032 == 1) && (slope032 == 0)) slope032 = 0.01;
    if ((do064 == 1) && (slope064 == 0)) slope064 = 0.01;
    if ((do128 == 1) && (slope128 == 0)) slope128 = 0.01;

    slopetotal = (slope001 + slope004 + slope016 + slope064) + (slope002 + slope008 + slope032 + slope128);  //diagonals have already been divided by 1.4 see above

    targets = do001 + do002 + do004 + do008 + do016 + do032 + do064 + do128;

    if (do002 == 1 || do008 == 1 || do032 == 1 || do128 == 1) thedistance = 1.41;

    // when pulling slope would be negative - need to make sure it is positive for proportions
    slope = ((double)selfhv - hv[client]) / thedistance; // placed back in Aug 30th 17

    //if ((slope<=0)|| (slope>6)) slope = 6;
    if ((slope == 0)) slope = 0.01;

    prop = 1.0 / targets; // default for SFD

 /* if ( targets > 1)
    {
    prop = (double) (slope/(slopetotal));
    } */

    return prop;
}
