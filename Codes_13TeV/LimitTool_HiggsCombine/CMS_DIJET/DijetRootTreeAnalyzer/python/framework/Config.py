import ConfigParser, os
import ROOT as rt

class Config(object):
    
    def __init__(self, fileName):
        if not os.path.exists(fileName):
            raise IOError("File not found: '%s'" % fileName)
        self.config = ConfigParser.ConfigParser()
        self.config.read(fileName)
    
    def __checkBox(self, box):
        if box not in self.config.sections():
            raise KeyError("The box '%s' was not found" % box)
        
    def get(self, box, var):
        result = None
        if self.config.has_option(box, var):
            result = self.config.get(box, var)
        elif self.config.defaults().has_key(var):
            result = self.config.defaults()[var]
        return result
    
    def has_option(self, box, var):
        result = False
        if self.config.has_option(box, var) or self.config.defaults().has_key(var):
            result = True
        return result
    
    def getVariables(self, box, lineTag='variables'):
        self.__checkBox(box)
        return eval(self.get(box,lineTag))
    
    def getVariablesRange(self, box, lineTag, workspace):
        #first define the variables
        workspace.defineSet(lineTag,'')  ###empty set with name corresponding to lineTag='variables'
        vars = self.getVariables(box, lineTag)  ###Variables in variables key in dijet.config file is in factory() input syntex, like x[5,-10,10]
        for v in vars:                          ###that is, var name with initial value and range.
            r = workspace.factory(v)
            workspace.extendSet(lineTag,r.GetName())  ##variables = ['mjj', 'th1x']

        #use a temporary RooWorkspace to process the ranges
        ws = rt.RooWorkspace('TMP')
        vars_ranges = self.getVariables(box, '%s_range' % lineTag)  ###['mjj_Low[453.,649.]', 'mjj_Blind[649.,838.]', 'mjj_High[838.,2037.]']
        for v in vars_ranges:
            a = ws.factory(v)
            
            #format must be VAR_RANGE
            name = a.GetName()  ##mjj_Low, mjj_Blind, mjj_High
            var_name, range_name = name.split('_') ###var_name = mjj, range_name = Low,Blind,High
            if workspace.var(var_name):
                workspace.var(var_name).setRange(range_name, a.getMin(), a.getMax())
            else:
                print "WARNING:: No variable found for range '%s'" % name
        return workspace
    
    def getRCuts(self, box):
        self.__checkBox(box)
        return eval(self.config.get(box,'rcuts'))
    
    def getBoxes(self):
        """Returns the names of the boxes defined in the config file"""
        if self.config.defaults().has_key('boxes'):
            return eval(self.config.defaults()['boxes'])
        return self.config.sections()
    
    def hasBinning(self, box):
        """Does the config for the box define a signal binning section?"""
        
        self.__checkBox(box)
        vars = self.getVariables(box)
        varNames = [v.split('[')[0] for v in vars]
        result = True
        for v in varNames:
            if not self.has_option(box, 'signal_%s' % v):
                result = False
                break
        return result
    
    def getBinning(self, box):
        """Returns the signal binning defined for the box"""
        self.__checkBox(box)
        vars = self.getVariables(box)  ###vars = ['mjj[453.,453.,2037.]','th1x[0,0,26]'], that is two variables. 
        
        varNames = [v.split('[')[0] for v in vars] ###runs over vars, and for each vars, return the first word before [, that is mjj and th1x.
        result = []
        for v in varNames:
            if self.has_option(box, 'signal_%s' % v):
                result.append(eval(self.get(box, 'signal_%s' %  v)) )
        return result        ###result contain binning for both mjj and th1x in two arrays within one array of two elements.
        
    def getPdfs(self, box, lineTag, workspace):
        pdfs = self.getVariables(box, lineTag)
        for p in pdfs:
            r = workspace.factory(p)
        return workspace
    
