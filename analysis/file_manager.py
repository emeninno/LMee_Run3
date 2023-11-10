import ROOT
from ROOT import TFile, TDirectory, TDirectoryFile, TList, THashList

class FileManager:
    def __init__(self,list_hierarchies):
        rootfile = TFile.Open(list_hierarchies[0],"READ"); 
        self.list = [];
        self.list.append(rootfile);
        ROOT.SetOwnership(rootfile,False);
        n = len(list_hierarchies);
        for i in range(1,n):
            #print(self.list[i-1].IsA());
            if self.list[i-1].IsA() == TFile.Class() or self.list[i-1].IsA() == TDirectory.Class() or self.list[i-1].IsA() == TDirectoryFile.Class():
                self.list.append(self.list[i-1].Get(list_hierarchies[i]));
            else:
                self.list.append(self.list[i-1].FindObject(list_hierarchies[i]));

        print(self.list);
        #print(self.list[-1].Print());

    def get(self,objname):
        if self.list[-1] is None:
            print("self.list[-1] is None. return None.");
            return None;
        else:
            if self.list[-1].IsA() == TFile.Class() or self.list[-1].IsA() == TDirectory.Class() or self.list[-1].IsA() == TDirectoryFile.Class():
                return self.list[-1].Get(objname);
            else:
                return self.list[-1].FindObject(objname);

    def close(self):
        if self.list[0].IsOpen():
            self.list[0].Close();
#_______________________________________________________________________
if __name__ == "__main__":

    list_hierarchies = [
    "AnalysisResults.root"
    ,"analysis-event-selection"
    ,"output"
    ];

    fm = FileManager(list_hierarchies);
    list_ev_after = fm.get("Event_AfterCuts");
    list_ev_after.ls();
    fm.close();
    #print(fm.get_list());
