import numpy as np
import math
import ROOT
from ROOT import TH1D, TH2D, TH3D, TMath, TH1F

def rebin_histogram_2d(h2org, nx,xmin,xmax,ny,ymin,ymax):
    #this supports only for constant bin width
    h2 = TH2D("{0}_rebin".format(h2org.GetName()),"{0}_rebin".format(h2org.GetName()), int(nx),xmin,xmax,int(ny),ymin,ymax); 

    h2.SetXTitle(h2org.GetXaxis().GetTitle());
    h2.SetYTitle(h2org.GetYaxis().GetTitle());

    bw_x_org = h2org.GetXaxis().GetBinWidth(1);
    bw_y_org = h2org.GetYaxis().GetBinWidth(1);
    bw_x = h2.GetXaxis().GetBinWidth(1);
    bw_y = h2.GetYaxis().GetBinWidth(1);
    if bw_x < bw_x_org :
        print("check bin width on x axis");
        return None;
    if bw_y < bw_y_org :
        print("check bin width on y axis");
        return None;
    if nx > h2org.GetNbinsX():
        print("check nbins x");
        return None;
    if ny > h2org.GetNbinsY():
        print("check nbins y");
        return None;
    if xmin < h2org.GetXaxis().GetXmin():
        print("check min. value of x");
        return None;
    if xmax > h2org.GetXaxis().GetXmax():
        print("check max. value of x");
        return None;
    if ymin < h2org.GetYaxis().GetXmin():
        print("check min. value of y");
        return None;
    if ymax > h2org.GetYaxis().GetXmax():
        print("check max. value of y");
        return None;

    #ng_x = int(bw_x / bw_x_org);
    #ng_y = int(bw_y / bw_y_org);
    #h2tmp = h2org.Clone("h2tmp");
    #h2tmp.RebinX(ng_x);
    #h2tmp.RebinY(ng_y);
    #now, bin width of h2org is adjusted for h2 on both x and y. !this is wrong 11.Feb.2022
    for ix in range(0,nx):
        for iy in range(0,ny):
            #h2.SetBinContent(ix+1, iy+1, h2tmp.GetBinContent(ix+1, iy+1));
            #h2.SetBinError(ix+1, iy+1, h2tmp.GetBinError(ix+1, iy+1));
            xlow  = h2.GetXaxis().GetBinLowEdge(ix+1);
            xhigh = h2.GetXaxis().GetBinLowEdge(ix+2);
            ylow  = h2.GetYaxis().GetBinLowEdge(iy+1);
            yhigh = h2.GetYaxis().GetBinLowEdge(iy+2);
            #print(xlow, xhigh, ylow, yhigh);
            bin_xlow_org  = h2org.GetXaxis().FindBin(xlow + 1e-6);
            bin_xhigh_org = h2org.GetXaxis().FindBin(xhigh - 1e-6);
            bin_ylow_org  = h2org.GetYaxis().FindBin(ylow + 1e-6);
            bin_yhigh_org = h2org.GetYaxis().FindBin(yhigh - 1e-6);
            n = 0;
            err = 0;
            #print(bin_xlow_org, bin_xhigh_org, bin_ylow_org, bin_yhigh_org);
            #print("=============");
            for ix2 in range(bin_xlow_org, bin_xhigh_org+1):
                for iy2 in range(bin_ylow_org, bin_yhigh_org+1):
                    #print(ix2, iy2);
                    n   += h2org.GetBinContent(ix2,iy2);
                    err += h2org.GetBinError(ix2,iy2) * h2org.GetBinError(ix2,iy2) ;
            #print(n,err);
            h2.SetBinContent(ix+1, iy+1, n);
            h2.SetBinError(ix+1, iy+1, TMath.Sqrt(err));
    return h2;
#______________________________________________________________________
def rebin_histogram(h1, arrX, isdiff, is_syst=False):
    h1tmp = h1.Clone("h1tmp");
    h1rebin = h1tmp.Rebin(len(arrX)-1,"h1rebin",arrX);
    h1rebin.SetName("{0}_rebin".format(h1.GetName()));
    h1rebin.Sumw2();

    if is_syst:
        for i in range(0,h1rebin.GetNbinsX()):
            x_min = h1rebin.GetBinLowEdge(i+1);
            x_max = h1rebin.GetBinLowEdge(i+2);
            bin0 = h1.FindBin(x_min + 1e-6);
            bin1 = h1.FindBin(x_max - 1e-6);
            y = 0;
            err = 0;
            for j in range(bin0, bin1+1):
                y += h1.GetBinContent(j);
                err += h1.GetBinError(j) / h1.GetBinContent(j) * h1.GetBinContent(j);
            h1rebin.SetBinContent(i+1, y);
            h1rebin.SetBinError(i+1, err);
            #print("check 0 rel syst. = ", h1rebin.GetBinError(i+1) / h1rebin.GetBinContent(i+1));

    if isdiff and (h1rebin.Class() == TH1D.Class() or h1rebin.Class() == TH1F.Class() ) : #do you want differential histogram? e.g. dN/dpT, dN/dm.
        h1rebin.Scale(1.,"width");
    return h1rebin;
#______________________________________________________________________
def slice_histogram(h2,x0,x1,axis,isdiff):
    h1 = 0;
    delta = 1e-6;
    if "x" in axis.lower():
        bin0 = h2.GetYaxis().FindBin(x0 + delta);
        bin1 = h2.GetYaxis().FindBin(x1 - delta);
        h1 = h2.ProjectionX("h1prjx_{0}".format(h2.GetName()),bin0,bin1,"");
    elif "y" in axis.lower():
        bin0 = h2.GetXaxis().FindBin(x0 + delta);
        bin1 = h2.GetXaxis().FindBin(x1 - delta);
        h1 = h2.ProjectionY("h1prjy_{0}".format(h2.GetName()),bin0,bin1,"");

    if isdiff and h1.Class() == TH1D.Class(): #do you want differential histogram? e.g. dN/dpT, dN/dm.
        h1.Scale(1.,"width");
    return h1;
#______________________________________________________________________
def slice_profile(h2,x0,x1,axis,isdiff=False):
    h1 = 0;
    delta = 1e-6;
    if "x" in axis.lower():
        bin0 = h2.GetYaxis().FindBin(x0 + delta);
        bin1 = h2.GetYaxis().FindBin(x1 - delta);
        h1 = h2.ProfileX("h1prfx_{0}".format(h2.GetName()),bin0,bin1,"");
    elif "y" in axis.lower():
        bin0 = h2.GetXaxis().FindBin(x0 + delta);
        bin1 = h2.GetXaxis().FindBin(x1 - delta);
        h1 = h2.ProfileY("h1prfy_{0}".format(h2.GetName()),bin0,bin1,"");

    if isdiff and h1.Class() == TProfile.Class(): #do you want differential profile?
        h1.Scale(1.,"width");
    return h1;
#______________________________________________________________________
def get_R_factor_2D(hULSnp_mix, hULSpn_mix, hLSpp_mix, hLSnn_mix):
    hR   = hULSnp_mix.Clone("hR");
    hR.Reset();
    hR.SetName("Rfactor_2D");
    hR.SetTitle("R factor in 2 dimensions");
    #hR.Sumw2();
    
    hULS_mix = hULSnp_mix.Clone("hULS_mix");#sum of np+pn
    if hULSpn_mix is not None:
        hULS_mix.Add(hULSpn_mix,1);
	
    Nm = hR.GetNbinsX();
    Pm = hR.GetNbinsY();
    #R = 1;
    #R_err = 0;
    for im in range(0,Nm):
        for ip in range(0,Pm):
            uls      = hULS_mix .GetBinContent(im+1, ip+1);
            uls_err  = hULS_mix .GetBinError(im+1, ip+1);
            lspp     = hLSpp_mix.GetBinContent(im+1, ip+1);
            lspp_err = hLSpp_mix.GetBinError(im+1, ip+1);
            lsnn     = hLSnn_mix.GetBinContent(im+1,ip+1);
            lsnn_err = hLSnn_mix.GetBinError(im+1,ip+1);
            R = 1.;
            R_err = 0.;
            if uls > 1e-6:
                if lspp * lsnn > 1e-6 :
                    R = uls / ( 2 * TMath.Sqrt( lspp*lsnn ) );
                    R_err = TMath.Sqrt( ( pow(lspp*lsnn_err*uls,2) + pow(lspp_err*lsnn*uls,2) + 4*pow(lspp*lsnn*uls_err,2) ) / ( 16 * pow(lspp*lsnn,3) ) );
                elif lspp + lsnn > 1e-6 :
                    R = uls / ( 2 * 0.5 * (lspp + lsnn) );
                    R_err = TMath.Sqrt( (pow(uls*lspp_err,2) + pow(uls*lsnn_err,2) + pow(lspp+lsnn,2)*pow(uls_err,2) ) / pow(lspp+lsnn,4) );
            else:
                R = 1.;
                R_err = 0.;
                #print("mix im+1 = {0} , x = {1} , R = {2} , uls = {3} , lspp = {4} , lsnn = {5}".format(im+1,h1R.GetBinCenter(im+1),R,uls,lspp,lsnn));
            hR.SetBinContent(im+1, ip+1, R);
            hR.SetBinError(im+1, ip+1, R_err);
    return hR;
#______________________________________________________________________
def get_R_factor(h1m_ULSnp_mix, h1m_ULSpn_mix, h1m_LSpp_mix, h1m_LSnn_mix):
    h1R   = h1m_ULSnp_mix.Clone("h1R");
    h1R.Reset();
    #h1R.Sumw2();
    h1R.SetName("Rfactor_1D");
    h1R.SetTitle("Rfactor in 1 dimension");

    h1m_ULS_mix = h1m_ULSnp_mix.Clone("h1m_ULS_mix");#sum of np+pn
    if h1m_ULSpn_mix is not None:
        h1m_ULS_mix.Add(h1m_ULSpn_mix,1.);

    Nm = h1R.GetNbinsX();
    #R = 1;
    #R_err = 0;
    for im in range(0,Nm):
        uls      = h1m_ULS_mix .GetBinContent(im+1);
        uls_err  = h1m_ULS_mix .GetBinError(im+1);
        lspp     = h1m_LSpp_mix.GetBinContent(im+1);
        lspp_err = h1m_LSpp_mix.GetBinError(im+1);
        lsnn     = h1m_LSnn_mix.GetBinContent(im+1);
        lsnn_err = h1m_LSnn_mix.GetBinError(im+1);
        R = 1.;
        R_err = 0.;

        if uls > 1e-6:
            if lspp * lsnn > 1e-6 :
                R = uls / ( 2 * TMath.Sqrt( lspp*lsnn ) );
                R_err = TMath.Sqrt( ( pow(lspp*lsnn_err*uls,2) + pow(lspp_err*lsnn*uls,2) + 4*pow(lspp*lsnn*uls_err,2) ) / ( 16 * pow(lspp*lsnn,3) ) );
            elif lspp + lsnn > 1e-6 :
                R = uls / ( 2 * 0.5 * (lspp + lsnn) );
                R_err = TMath.Sqrt( (pow(uls*lspp_err,2) + pow(uls*lsnn_err,2) + pow(lspp+lsnn,2)*pow(uls_err,2) ) / pow(lspp+lsnn,4) );
        else:
            R = 1.;
            R_err = 0.;
        #print("mix im+1 = {0} , x = {1} , R = {2} , uls = {3} , lspp = {4} , lsnn = {5}".format(im+1,h1R.GetBinCenter(im+1),R,uls,lspp,lsnn));
        h1R.SetBinContent(im+1, R);
        h1R.SetBinError(im+1, R_err);
    return h1R;

#______________________________________________________________________
def get_bkg1D(hLSpp, hLSnn):
    hbkg = hLSpp.Clone("hbkg");
    hbkg.Reset();
    nbins = hbkg.GetNbinsX();
    for im in range(0,nbins):
        lspp = hLSpp.GetBinContent(im+1);
        lspp_err = hLSpp.GetBinError(im+1);
        lsnn = hLSnn.GetBinContent(im+1);
        lsnn_err = hLSnn.GetBinError(im+1);

        if lspp * lsnn > 1e-6 : # use geometric mean....need to check this
            hbkg.SetBinContent(im+1,2.0 * TMath.Sqrt(lspp * lsnn));
            hbkg.SetBinError(im+1,2.0 * TMath.Sqrt(lspp_err + lsnn_err));
        elif lspp + lsnn > 1e-6 :# use arithmetic mean to avoid empty bin issus
            hbkg.SetBinContent(im+1,2.0 * (lspp * lsnn));
            hbkg.SetBinError(im+1,2.0 * TMath.Sqrt(lspp_err + lsnn_err));
	   
    return hbkg;
												
#______________________________________________________________________
def get_bkg2D(hLSpp, hLSnn):
    hbkg = hLSpp.Clone("hbkg");
    hbkg.Reset();
    nBinsX = hbkg.GetNbinsX();
    nBinsY = hbkg.GetNbinsY();

    for im in range(0,nBinsX):
        for ip in range(0,nBinsY):
            lspp = hLSpp.GetBinContent(im+1,ip+1);
            lspp_err = hLSpp.GetBinError(im+1,ip+1);
            lsnn = hLSnn.GetBinContent(im+1,ip+1);
            lsnn_err = hLSnn.GetBinError(im+1,ip+1);

            if lspp * lsnn > 1e-6 : # use geometric mean....need to check this
               	hbkg.SetBinContent(im+1,ip+1,2.0 * TMath.Sqrt(lspp * lsnn));
               	hbkg.SetBinError(im+1,ip+1,2.0 * TMath.Sqrt(lspp_err + lsnn_err));
            elif lspp + lsnn > 1e-6 : # use arithmetic mean to avoid empty bin issus  		
                hbkg.SetBinContent(im+1,ip+1,2.0 * (lspp * lsnn));
                hbkg.SetBinError(im+1,ip+1,2.0 * TMath.Sqrt(lspp_err + lsnn_err));
            
    return hbkg;
    	
#______________________________________________________________________
def get_corrected_bkg(h1R, h1m_LSpp_same, h1m_LSnn_same):
    h1bkg = h1m_LSpp_same.Clone("h1bkg");
    h1bkg.Reset();
    #h1bkg.Sumw2();

    Nm = h1bkg.GetNbinsX();
    for im in range(0,Nm):
        lspp     = h1m_LSpp_same.GetBinContent(im+1);
        lspp_err = h1m_LSpp_same.GetBinError(im+1);
        lsnn     = h1m_LSnn_same.GetBinContent(im+1);
        lsnn_err = h1m_LSnn_same.GetBinError(im+1);
        R        = h1R.GetBinContent(im+1);
        R_err    = h1R.GetBinError(im+1);

        #print("im+1 = {0} , x = {1} , R = {2} , lspp = {3} , lsnn = {4}".format(im+1,h1R.GetBinCenter(im+1),R,lspp,lsnn));

        bkg = 0;
        bkg_err = 0;
        if lspp * lsnn > 1e-6 : #geometric mean
            bkg = 2.0 * R * TMath.Sqrt(lspp * lsnn);
            bkg_err = TMath.Sqrt( pow(R_err/R,2) + 1./4*pow(lspp_err/lspp,2) + 1./4*pow(lsnn_err/lsnn,2) ) * bkg;
            #print("im+1 = {0} , x = {1} , R = {2} , lspp = {3} , lsnn = {4}, bkg_err = {5}".format(im+1,h1R.GetBinCenter(im+1),R,lspp,lsnn,bkg_err));
            #bkg_err = TMath.Sqrt( R*R * ( pow(lspp*lsnn_err,2) + pow(lsnn*lspp_err,2) ) / (lspp*lsnn) );
            #print("old im+1 = {0} , x = {1} , R = {2} , lspp = {3} , lsnn = {4}, bkg_err = {5}".format(im+1,h1R.GetBinCenter(im+1),R,lspp,lsnn,bkg_err));
        elif lspp + lsnn > 1e-6 : #arithmetic mean
            bkg = 2.0 * R * 0.5 * (lspp + lsnn);
            bkg_err = TMath.Sqrt( pow(R_err/R,2) + (pow(lspp_err,2)+pow(lsnn_err,2))/pow(lspp+lsnn,2) ) * bkg;
            #print("im+1 = {0} , x = {1} , R = {2} , lspp = {3} , lsnn = {4}, bkg_err = {5}".format(im+1,h1R.GetBinCenter(im+1),R,lspp,lsnn,bkg_err));
            #bkg_err = TMath.Sqrt( R*R * (lspp_err*lspp_err + lsnn_err*lsnn_err) );
            #print("old im+1 = {0} , x = {1} , R = {2} , lspp = {3} , lsnn = {4}, bkg_err = {5}".format(im+1,h1R.GetBinCenter(im+1),R,lspp,lsnn,bkg_err));
        h1bkg.SetBinContent(im+1,bkg);
        h1bkg.SetBinError(im+1,bkg_err);

    return h1bkg;
#______________________________________________________________________
def get_corrected_bkg_simple(R, R_err, h1m_LSpp_same, h1m_LSnn_same):
    h1bkg = h1m_LSpp_same.Clone("h1bkg");
    h1bkg.Reset();
    #h1bkg.Sumw2();

    Nm = h1bkg.GetNbinsX();
    for im in range(0,Nm):
        lspp     = h1m_LSpp_same.GetBinContent(im+1);
        lspp_err = h1m_LSpp_same.GetBinError(im+1);
        lsnn     = h1m_LSnn_same.GetBinContent(im+1);
        lsnn_err = h1m_LSnn_same.GetBinError(im+1);

        bkg = 0;
        bkg_err = 0;
        if lspp * lsnn > 1e-6 :
            bkg = 2.0 * R * TMath.Sqrt(lspp * lsnn);
            bkg_err = TMath.Sqrt( R*R * ( pow(lspp*lsnn_err,2) + pow(lsnn*lspp_err,2) ) / (lspp*lsnn) );
        elif lspp + lsnn > 1e-6 : #arithmetic mean
            bkg = 2.0 * R * 0.5 * (lspp + lsnn);
            bkg_err = TMath.Sqrt( R*R * (lspp_err*lspp_err + lsnn_err*lsnn_err) );
        h1bkg.SetBinContent(im+1,bkg);
        h1bkg.SetBinError(im+1,bkg_err);

    return h1bkg;
#______________________________________________________________________
def get_bkg_subtracted(h1m_ULS_same, h1bkg):
    h1sig = h1m_ULS_same.Clone("h1sig");
    #h1sig.Sumw2();
    h1sig.Add(h1bkg,-1);
    return h1sig;
#______________________________________________________________________
def get_SBratio(h1sig,h1bkg):
    h1r = h1sig.Clone("h1r");
    h1r.Reset();
    h1r.Divide(h1sig,h1bkg,1.,1.,"B");
    return h1r;

#______________________________________________________________________
#______________________________________________________________________
def get_significance(h1sig,h1bkg):
    h1r = h1sig.Clone("h1r");
    h1r.Reset();

    n = h1r.GetNbinsX();
    for i in range(0,n):
        s = h1sig.GetBinContent(i+1);
        b = h1bkg.GetBinContent(i+1);
        s_err = h1sig.GetBinError(i+1);
        b_err = h1bkg.GetBinError(i+1);
        sig = s/sqrt(s + 2.*b);
        sig_err = sqrt( pow(pow(s+2*b,-1/2.) - s/2.*pow(s+2.*b, -3/2) ,2) *pow(s_err,2) + pow( 2*s * pow(s+2*b,-3/2.) ,2) * pow(b_err,2)  );
        h1r.SetBinContent(i+1,sig);
        h1r.SetBinError(i+1,sig_err);
    return h1r;

#______________________________________________________________________
#______________________________________________________________________
#______________________________________________________________________
#______________________________________________________________________
