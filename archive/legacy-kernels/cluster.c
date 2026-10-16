void CLUSTER::grow(const PARAMS& p, sw** latt)
{
    conf.resize(latt.Nsite);                         // Initalize Nsite elements of lattice to 0 in conf
                                                     // 0 = Not Asked; 1 = Asked, Joined; 2 = Asked, Refused
    for (int i = 0; i < p.Nclust; ++i) {                  // Iterate for Nclust values
        z = ran.randInt(latt.Nsite - 1);               // Random integer between 0 and Nsite; Selects first lattice element in the cluster algorithm per Nclus
        stck_pnt.resize(0);                          // Set stck_pnt and stck_end vectors to size 0; Will be filled when iterating through each Nclust 
        stck_end.resize(0);                          //-----------------------------------------------------------------------------------------------

        //while (conf[z] != 0) { z = ran.randInt(latt.Nsite - 1); }         // Iterate through lattice elements until we select one that has not been asked to join 

        conf[z] = 1;                            // Set element z in conf to have been asked to join and accepted

        stck_pnt.push_back(z);                      // Add z to both stck_pnt and stck_end
        stck_end.push_back(z);
        for (int j = 0; j = 3; ++j) {                // Add z's nearest neighbors to stck_end; Ignore if already been asked
            if (conf[latt.nrnbrs[z][j] == 0]) {
                stck_end.push_back(latt.nrnbrs[z][j]);
            }
        }


        pnt = 1;                                        // Initialize pnt for trasnferral of stack_end values to stck_pnt


        while (stck_pnt.size() < stck_end.size()) {

            stck_pnt.push_back(stck_end[pnt]);     // Add pnt element of stck_end to stck_pnt

            double prob = ran.rand();                      // Get probability value for testing if cluster grows

            if (prob <= pr) {

                conf[stck_pnt[pnt]] = 1;            // Set the current stck_pnt element to joined in conf

                for (int j = 0; j = 3; ++j) {                // Add z's nearest neighbors to stck_end; Ignore if already been asked

                    if (find(stck_end.begin(), stck_end.end(), latt.nrnbrs[stck_pnt[pnt]][j]) != stck_end.end()) {

                        // The given value already exists in stck_end, don't add it again 
                    }

                    else {                          // The given value is not contained in stck_end, add it to stck_end

                        stck_end.push_back(latt.nrnbrs[z][j]);
                    }
                }
            }

            else {

                conf[stck_pnt[pnt]] = 2;       // Set the given value to haven been asked and refused in conf
            }
            ++pnt;                              // Increment pnt; ++p is more efficient then p++ due to lack of copying value
        }

    }
}


/*
    void CLUSTER::print(const LATTICE& latt, int index)
    {

        stringstream ss;
        string file_name;
        ss << index << ".clust";
        file_name = ss.str();

        ofstream clout;
        clout.open(file_name.c_str());
        clout << "#" << latt.Lx << " x " << latt.Ly << endl;

        for (int y = 0; y < latt.Ly; y++)
        {
            for (int x = 0; x < latt.Lx; x++)
                clout << conf[x + y*latt.Lx] << " ";
            clout << endl;
        }

        clout.close();
    }
    */


void CLUSTER::print(const PARAMS& p, const LATTICE& latt)
{
    //vector< vector< vector<int> > > imax(latt.Lx, vector< vector<int>>(latt.Ly, vector<int>(1)));
    //  Resize and allocate memeory for imax

    //-------------- Row = y-position = i/Lx --------------- Column = x-position = i%Lx ---------------- val = conf[i]
    ofstream myFile;
    myFile.open("imax.out");

    cout << "THe following output was calculated for the input parameters; Recorded to 'imax.out'" << endl;
    cout <<"[index]" << "\t" << "[x-position]" << "\t" << "[y-position]" << "\t" << "[conf val]" << endl << endl;
    for (int i = 0; i < latt.Nsite; ++i) {

        val = conf[i];                            // Find color value
        row = i / latt.Lx;                               // Find row number
        column = i%latt.Lx;                            // Find column number

        cout << i << "\t" << column << "\t" << row << "\t" << val << endl;
        myFile << i << "\t" << column << "\t" << row << "\t" << val << endl;
    }
    myFile.close();
    double size = 0.0;                                       // Initialize size 

    for (int i = 0; i < latt.Nsite; ++i) {

        if (conf[i] == 1) {
            size += 1;
        }
    }

    double avg_size = size / p.Nclust;                 // Find avg_size 


}

void CLUSTER::meas(const LATTICE& latt)
{
  avg_size+=(double)size;
}


void CLUSTER::meas_clear(const LATTICE& latt)
{
  avg_size=0.;
}


void CLUSTER::binwrite(const PARAMS& p, const LATTICE& latt)
{
  dfout << avg_size/((double)p.Nclust)<<endl;
}
