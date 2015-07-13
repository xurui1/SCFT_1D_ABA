string IntToStr(int n)
{
    stringstream result;
    result << n;
    return result.str();
}

void FreeEnergy(double **w, double **phi, double *eta, int *Ns, double ds, double *chi, double dr, double **chiMatrix, double *mu, double fE_hom, double volume, double *f){
    
    
    double  currentfE, oldfE, deltafE;
    int     maxIter=100000;
    double precision=1e-5;          //convergence condition
    int     i,x,iter,ii,jj,radius;
    int     mmb;
    double  Q;
    double  fE_int, fES;            //interaction free energy and chain partition function fE
    double  epsilon, gamma;
    double  *delphi;
    double *sigma;
    double  **delW;
    double  **newW;
    double  deltaW;
    
    //Arrays for updating the omega fields
    delW=create_2d_double_array(ChainType,Nr,"delW");
    delphi=create_1d_double_array(Nr,"delphi");
    sigma=create_1d_double_array(Nr,"sigma");
    newW=create_2d_double_array(ChainType,Nr,"newW");
    
    currentfE=0.0;
    deltafE=0.0;
    
    epsilon=0.05;
    gamma=0.05;
    
    iter=0;
    mmb=1;
    std::ofstream outputFile1("./results/fE.dat");
    std::ofstream outputFile2("./results/fE_R.dat");
    std::ofstream outputFile3("./results/fE_Chi.dat");
    
    /*for (x=0;x<25;x++){*/
    
    /*for (radius=0;radius<150;radius++){
        volume=vol(dr);
        if (radius>1){omega(w);}*/
    for (iter=0;iter<maxIter;iter++){
        
        fE_int=0.0;
        fES=0.0;
        deltaW=0.0;

        
        Q=Conc(phi,w,Ns,ds,dr,mu,volume);          //Calculate Chain partition function for both AB and C
        
        
        Incomp(eta,phi,delphi);              //Enforce incompressibility condition
        output(dr,phi,w);                   //Output some data to file
        
        if (mmb==1){
            Pin(sigma, phi);
        }
        

        
        //Calculate components for new field and interaction free energies
        for(i=0;i<Nr;i++){
            for(ii=0;ii<ChainType;ii++){
                newW[ii][i]=0.0;            //set field update to zero
                for(jj=0;jj<ChainType;jj++){
                    newW[ii][i]+=(chiMatrix[ii][jj]*phi[jj][i]);
                }
                newW[ii][i]+=eta[i];
                
                if (mmb==1){
                    if (ii==0 || ii == 2 || ii==4){
                        newW[ii][i]-=sigma[i];
                    }
                    else if (ii==1 || ii==3){
                        newW[ii][i]+=sigma[i];
                    }
                }
                delW[ii][i]=newW[ii][i]-w[ii][i];
                w[ii][i]+=(gamma*delW[ii][i]-epsilon*delphi[i]);     //update omega field
                deltaW+=fabs(delW[ii][i])*dV(i,dr);
                }
        }
        fE_int=fE(newW,phi,chiMatrix,dr,volume);
        
        //Normalize by box size
        deltaW/=volume;
        
       
        //Update free energy
        fES=Q;
        oldfE=currentfE;
        currentfE=-fES+fE_int;
        
       /* if (fabs(currentfE)>1e3){
            cout<<"fE too large"<<endl;
            exit(EXIT_FAILURE);
        }*/
        deltafE=fabs(currentfE-oldfE);
        
        //Print free energy, difference in free energy, change in omega field to screen
        std::cout<<iter<<" fE:"<<currentfE<< " dfE:"<<currentfE-fE_hom<<" " << deltaW<<" "<<fE_hom<<std::endl;
        outputFile1 << iter << " " << currentfE<< " " << currentfE-fE_hom<<" "<< deltaW<<std::endl;
        

        if (deltafE<precision && deltaW<precision){break;} //Convergence condition
        
    }
        
    /*************************Loop for determining fE as f(Chi)***************************************/
    /*    outputFile3 << x <<" "<<chi[0]<<" "<< chi[1]<<" "<<chi[2]<<" "<<currentfE-fE_hom<<std::endl;
        chi[0]-=1.0;
        chi[1]-=1.0;
        Xmatrix(chiMatrix,chi);
        fE_hom=homogfE(mu,chiMatrix,f);

    }*/
    /**************************************************************************************************/
    
    
    /**************************Loop for determining fE as f(radius)***********************************/
     /*  outputFile2 << r_0 << " "<<currentfE-fE_hom<<std::endl;
        //create name
        ofstream outFile;
        string filename;
        filename="./results/phi_r_" + IntToStr((int)r_0) + ".dat"; // C++11 for std::to_string
        //create file
        outFile.open(filename.c_str());
        
        for (i=0;i<Nr;i++){
            outFile <<i*dr<<" "<<phi[0][i]<<" "<<phi[1][i]<<" "<<phi[2][i]<<" "<<phi[3][i]<<" "<<phi[4][i]<<" "<<phi[5][i]<<std::endl;;
        }
        outFile.close();
        
        if (r_0<20){
            r_0+=1.0;
        }
        else if (r_0<50){
            r_0+=2.0;
        }
        else if (r_0<100){
            r_0+=5.0;
        }
        else{
            r_0+=10.0;
        }
        
    }*/
    /**************************************************************************************************/
      
      
    outputFile1.close();
    outputFile2.close();
    outputFile3.close();
    
    
    destroy_1d_double_array(delphi);
    destroy_2d_double_array(delW);
    destroy_2d_double_array(newW);
    
}
