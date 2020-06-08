#!/bin/env python
#
"""
 Interface for configuring pre/post icetray scripts

 copyright  (c) 2005 the icecube collaboration

 @version: $Revision: $
 @date: $Date: $
 @author: Juan Carlos Diaz Velez <juancarlos@icecube.wisc.edu>
"""

import os
import re
import sys
import math
import time
import string
import shutil
import logging


class IPBaseClass:
    """
    This class provides an interface for preprocessing files in iceprod
    """

    def __init__(self):
        
        self.description   = {}
        self.AddParameter("execute","boolean condition to execute", True)

        self.status        = 0
        self.configured    = False
        self.parser        = None

        # Aggregate CPU times
        self.realtime      = 0.0
        self.usertime      = 0.0
        self.systime       = 0.0

        self.logger = logging.getLogger(self.__class__.__name__)
    
    @property
    def parameters(self):
        return self.__dict__
    
    def Configure(self,tray):
        self.logger.info("Configuring IceProd module '%s'" % self.__class__.__name__)
        self.configured    = True

    def SetParser(self, parser):
        """
        Set the ExpParser object
        """
        self.parser        = parser

    def AddParameter(self, param,description, value,type=None):
        """
        Add parameter value for plugin
        """
        setattr(self,param.lower(),value)
        self.description[param.lower()] = description

    def GetParameter(self, param):
        """
        Get parameter value for plugin
        """
        if param.lower() not in self.description:
            raise Exception("Attemting to get parameter %s not added by %s" % (param,self))
        return getattr(self,param.lower())

    def SetParameter(self, param, value):
        """
        Set parameter value for plugin
        """
        if param.lower() not in self.description:
            self.logger.error(self.ShowParameters())
            raise Exception("param %s was configured but not added by %s" % (param,self))

        setattr(self,param.lower(),value)
        self.logger.info("%s:%s" %(param,value))
        return self.SetParameter

    def SetParameters(self, *args):
        """
        Set multiple parameters for plugin
        """
        for key,value in args:
            self.SetParameter(key,value)
        return self


    def Execute(self,stats):
        self.logger.info("execute %s: %s" % (self.__class__.__name__,self.GetParameter("execute")))
        return self.GetParameter("execute")

    def ShowParameters(self):
        return zip(
            self.parameters.keys(),
            map(self.__dict__.get,self.__dict__.keys()),
            map(self.description.get,self.parameters.keys())
            )

    def Finish(self,stats={}):
        self.logger.info("finish %s: %s" % (self.__class__.__name__,self.GetParameter("execute")))
        return 0

    def PrintStats(self,stats={}):
        header = "Stats produced by " + str(self.__class__) +":"
        self.logger.info(header)
        self.logger.info("-"*len(header))
        for k,v in stats.items():
            self.logger.info("%s,%s" %(k,v))

class Hello(IPBaseClass):

    def __init__(self):
        IPBaseClass.__init__(self)
        self.AddParameter("greeting","String to write to screen","Hello World")

    def Execute(self,stats):
        if not IPBaseClass.Execute(self,stats):
            return 0
        self.logger.info(self.GetParameter("greeting"))
        return 0

class GoodBye(IPBaseClass):

    def __init__(self):
        IPBaseClass.__init__(self)
        self.AddParameter("greeting","String to write to screen","Good bye World")

    def Execute(self,stats):
        if not IPBaseClass.Execute(self,stats):
            return 0
        self.logger.info(self.GetParameter("greeting"))
        return 0

class GenericModule(IPBaseClass):
    def __init__(self):
        IPBaseClass.__init__(self)
        self.AddParameter("Executable","Name of executable to run","")
        self.AddParameter("Arguments","string of arguments to pass to executable",[])

    def Execute(self,stats):
        if not IPBaseClass.Execute(self,stats): return 0
        executable = self.GetParameter("Executable")
        arguments  = self.GetParameter("Arguments")
        cmd = "%s %s" % (executable," ".join(arguments))
        self.logger.info("executing:", cmd)
        return os.system(cmd)

def string_splitter_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

def handle_bool(value):
  if value: return "store_true"
  return "store_false"
        
try:
   from optparse import OptionParser
except:
   logging.error("could not import option parser. skipping definitions of class ParsingModule")
else:
   import types

   class ParsingModule(IPBaseClass):
        def __init__(self,opt_parser=None):
           if not opt_parser:
              opt_parser = OptionParser()
           self._opt_parser = opt_parser
           IPBaseClass.__init__(self)
           

        def AddParameter(self, param, description, value, explicit_type=None):
           try:
              IPBaseClass.GetParameter(self,param)
           except Exception:
              pass
           else:
              raise Exception('Tried to add parameter %s twice'%param)
           IPBaseClass.AddParameter(self,param,description,value)
           if explicit_type is not None:
              if explicit_type == 'list':
                  self._opt_parser.add_option('', '--%s' % param, dest=param.lower(), help=description, 
                      default=value, type='string', action='callback', callback=string_splitter_callback)
              else:
                  self._opt_parser.add_option('', '--%s' % param, dest=param.lower(),
                              type=explicit_type, help=description, default=value )
           elif type(value) == int:
              self._opt_parser.add_option('', '--%s' % param, dest=param.lower(),
                              type="int", help=description, default=value )
           elif type(value) == float:
              self._opt_parser.add_option('', '--%s' % param, dest=param.lower(),
                              type="float", help=description, default=value )
           elif type(value) == bool:
              if value is True:
                  self._opt_parser.add_option('', '--no-%s' % param, dest=param.lower(), 
                             action='store_false', help=description, default=value )
              else:
                  self._opt_parser.add_option('', '--%s' % param, dest=param.lower(), 
                             action='store_true', help=description, default=value )
           elif type(value) == list:
              self._opt_parser.add_option('', '--%s' % param, dest=param.lower(), help=description, 
                      default=value, type='string', action='callback', callback=string_splitter_callback)
           else:
              self._opt_parser.add_option('', '--%s' % param, type='string', dest=param.lower(), help=description, default=value )

        def Execute(self, stats={}):
            if not IPBaseClass.Execute(self,stats): return False
            return True

        def ExecuteOpts(self, stats={},opts = None):
            if opts is None: opts = self._opt_parser.parse_args()
            (self.parser_opts, self.parser_args) = opts
            for key in self.description.keys():
	            self.SetParameter(key,getattr(self.parser_opts,key.lower()))
            retval = self.Execute(stats)
            self.PrintStats(stats)
            return retval


   class FromOptionParser(IPBaseClass):
        """Converts an OptionParser object to an IceProd module.
           Takes a parser and a main function as input.
           Will pass a dict with the name,value pairs to the main function."""
        def __init__(self,parser,main):
            IPBaseClass.__init__(self)
            self.main = main
            self.args = parser.defaults
            for name in self.args:
                value = self.args[name]
                description = ''
                IPBaseClass.AddParameter(self,name,description,value)
        
        def Execute(self, stats={}):
            if not IPBaseClass.Execute(self,stats): return 0
            options = {}
            for name in self.args:
                options[name] = self.GetParameter(name)
            ret = self.main(options,stats)
            if ret is None:
                return 0
            else:
                return ret
