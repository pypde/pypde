import sys
from pyparsing \
    import nums, alphas, alphanums, quotedString, cStyleComment,cppStyleComment,           restOfLine, \
           delimitedList, removeQuotes, CaselessLiteral, Word, Dict, \
           Optional, ZeroOrMore, OneOrMore, Combine, Group, Suppress, Forward \
           , nestedExpr, dictOf
import mesh

##------------------------------------------------------------------------------
class CaseManager(object):
    def __init__(self,prefix,filename):

        filename = ''.join([prefix,'/',filename]);
        caseDict = self.__parseInput(filename)        

        self.debugDict = {}
        self.probesDict = {}

        self.problemDict = caseDict.get('problem')
        self.machineDict = caseDict['machine']
        self.modelDict = caseDict['problem']['model']
        self.timeDict = caseDict['problem']['discretization']['time']
        self.spaceDict = caseDict['problem']['discretization']['space']
        self.outputDict = caseDict['output']
        self.condsDict = caseDict['problem']['conditions']
        if 'debug' in caseDict:
            self.debugDict = caseDict['debug']
        if 'probes' in caseDict:
            self.probesDict = caseDict['probes']

        #self.prefix = self.machineDict.get('prefix')
        self.prefix = prefix
        fDomain = ''.join([prefix,'/',self.spaceDict.get('domainFile')])
        fBC = ''.join([prefix,'/',
                       self.condsDict.get('boundary').get('bcFile')])

        # get boundary conditions to build bc's later
        self.bcDict = self.__parseBCs(fBC)

        # get mesh
        self.msh = mesh.Mesh(fDomain)

    ##--------------------------------------------------------------------------

    ##--------------------------------------------------------------------------
    def __parseBCs(self,filenameBC):
        fB = open(filenameBC,'r')
        expressions = fB.read();

        EQ,LBRACE,RBRACE,SEMI = map(Suppress,"={};")
        LPAREN,RPAREN = map(Suppress,'()')

        intEntry = Combine(Optional('-') + Word(nums)) \
                   .setParseAction(lambda t : int(t[0]))
        realEntry = Combine(Optional('-') + Optional(Word(nums)) \
                            + '.' + (Word(nums))) \
                   .setParseAction(lambda t : float(t[0]))
        realEntry2 = Combine(Optional('-') + (Word(nums)) \
                            + '.' + Optional(Word(nums))) \
                   .setParseAction(lambda t : float(t[0]))
        numEntry = realEntry | realEntry2 | intEntry
        sciEntry = Combine(  numEntry \
                           + (CaselessLiteral('E') \
                              + Word('+-'+nums,nums))) \
                 .setParseAction(lambda t : float(t[0]))
        realOrSciEntry = sciEntry | realEntry | realEntry2

        # define tokens, expressions and entries
        bcObjEntry = Forward()
        bcObjListEntry = Forward()
        #listEntry = Forward()

        keyToken = Word(alphas+"_", alphanums+"_")
        entryToken = ( keyToken | realOrSciEntry | intEntry | 
                       quotedString.copy().setParseAction(removeQuotes))

        # define lists
        #numList = Group(OneOrMore(realOrSciEntry))
        #intList = Group(OneOrMore(intEntry))
        #keyTokenList = Group(OneOrMore(keyToken))
        numList = (OneOrMore(realOrSciEntry)).setParseAction(lambda t : list([t[:]]))
        intList = (OneOrMore(intEntry)).setParseAction(lambda t : list([t[:]]))
        keyTokenList = (OneOrMore(keyToken)).setParseAction(lambda t : list([t[:]]))
        entryList = numList | intList | keyTokenList

        bcObjExpr = Group(keyToken + bcObjEntry) \
                | Group (keyToken + LBRACE + RBRACE)
        bcObjListExpr = Group(keyToken + bcObjListEntry)
        #listExpr = (  Group(keyToken + listEntry) \
        #            | Group(keyToken + LPAREN + RPAREN))
        expr = Group(keyToken + entryToken + SEMI) \
             | Group(keyToken + LPAREN + entryList + RPAREN + SEMI)
            # | listExpr + SEMI

        # expr needs to come before the list ones
        #bcMixedExpr = (bcObjExpr | expr | listExpr)
        bcMixedExpr = (bcObjExpr | expr | bcObjListExpr)

        bcObjListEntry << ( LPAREN + Dict(OneOrMore(bcObjExpr)) + RPAREN)
        #listEntry <<   ( LPAREN + entryList + RPAREN)
        #             | ( LPAREN + OneOrMore(listExpr) + RPAREN))#  this last line allows nesting 
        bcObjEntry << (  LBRACE + Dict(OneOrMore(bcMixedExpr)) + RBRACE)

        #print expressions
        bcDict_t = Dict(OneOrMore(bcMixedExpr)) \
                  .ignore(cStyleComment).ignore(cppStyleComment)\
                  .parseString(expressions,parseAll=True)

        fB.close()

        #print bcDict_t.dump()
        #bcTypeList = []
        #for i in range(len(msh.patchList)):
        #    #print bcDict_t[msh.patchList[i][0]]
        #    bcTypeList.append(bcDict_t[msh.patchList[i][0]].keys())
 
        #bcDict = {}
        #for iB in bcDict_t.keys():
        #    subBCDict = {}
        #    for kkk in range(len(bcDict_t[iB])):
        #        print kkk
        #        #print bcDict_t[iB][kkk][1]
        #        #print bcDict_t[iB][kkk][1][1]
        #        print bcDict_t[iB].values()
        #        subBCDict.update({bcDict_t[iB][kkk][0]:bcDict_t[iB][kkk][1]})
        #    bcDict[iB] = subBCDict.copy()
        #    #bcDict[iB] = subBCDict.deepcopy()

        #print bcDict
        #sys.exit(0)
        return bcDict_t
    ##--------------------------------------------------------------------------

    ##--------------------------------------------------------------------------
    def __parseInput(self,filename):

        ##----------------------------------------------------------------------
        ## Read input file
        ##----------------------------------------------------------------------
        print " ... Parsing input file \"", filename ,"\""
        sys.stdout.flush()

        f = open(filename,'r')
        expressions = f.read()
        #print expressions
        
        EQ,LBRACE,RBRACE,SEMI = map(Suppress,'={};')
        LPAREN,RPAREN = map(Suppress,'()')
        
        intEntry = Combine(Optional('-') + Word(nums)) \
                   .setParseAction(lambda t : int(t[0]))
        realEntry = Combine(Optional('-') + Optional(Word(nums)) \
                            + '.' + (Word(nums))) \
                   .setParseAction(lambda t : float(t[0]))
        realEntry2 = Combine(Optional('-') + (Word(nums)) \
                            + '.' + Optional(Word(nums))) \
                   .setParseAction(lambda t : float(t[0]))
        numEntry = realEntry | realEntry2 | intEntry

        sciEntry = Combine(  numEntry \
                           + (CaselessLiteral('E') \
                              + Word('+-'+nums,nums))) \
                 .setParseAction(lambda t : float(t[0]))
        
        realOrSciEntry = sciEntry | realEntry | realEntry2

        # define tokens, expressions and entries
        objEntry = Forward()
        listEntry = Forward()
        objListEntry = Forward()
        flListEntry = Forward()

        keyToken = Word(alphas+'_', alphanums+'_')
        entryToken = ( keyToken | realOrSciEntry | intEntry | 
                       quotedString.copy().setParseAction(removeQuotes))

        # define lists
        numList = OneOrMore(realOrSciEntry).setParseAction(lambda t : list([t[:]]))
        intList = OneOrMore(intEntry).setParseAction(lambda t : list([t[:]]))
        keyTokenList = OneOrMore(keyToken).setParseAction(lambda t : list([t[:]]))
        entryList = numList | intList | keyTokenList

        objExpr = (  Group(keyToken + objEntry)
                   | Group(keyToken + LBRACE + RBRACE))
        objListExpr = Group(keyToken + objListEntry)
        listExpr = (  Group(keyToken + listEntry) \
                    | Group(keyToken + LPAREN + RPAREN))
        expr = Group(keyToken + entryToken + SEMI) \
             | listExpr + SEMI

        # expr needs to come before the list ones
        mixedExpr = (objExpr | expr | objListExpr | listExpr)
        
        objListEntry << ( LPAREN + Dict(OneOrMore(objExpr)) + RPAREN)
        listEntry <<   ( LPAREN + entryList + RPAREN)
                      #| ( LPAREN + OneOrMore(listExpr) + RPAREN))#  this last line allows nesting 


        objEntry << Group(LBRACE + Dict(OneOrMore(mixedExpr)) + RBRACE)
        
        caseDict = Dict(OneOrMore(mixedExpr)) \
                   .ignore(cStyleComment).ignore('//' + restOfLine) \
                   .parseString(expressions,parseAll=True)

        print '     Done'
        sys.stdout.flush()
        #print ' ... Input file dictionary:'
        #print caseDict.dump()
        #sys.stdout.flush()
        f.close()


        ##----------------------------------------------------------------------
        ## shorter parsing code (too flexible)
        ##----------------------------------------------------------------------

        #print 'NEW PARSING'

        #f = open('zpyfInput','r')
        #expressions = f.read()
        #value = Forward()
        #result = Forward()

        #testDict = nestedExpr('{','}',result)
        ##testDict = 
        ##    Group (LBRACE + Optional(Dict(delimitedList(result,delim=SEMI))) 
        ##           + Optional(SEMI) + RBRACE)
        #testList = Group (LPAREN + Optional(delimitedList(result,delim=SEMI)) 
        #                  + Optional(SEMI) + RPAREN)
        #testList2 = Group (LPAREN + Optional(OneOrMore(entryList)) 
        #                   + Optional(SEMI) + RPAREN)

        #result << Group ((keyToken | intEntry) + value )
        #value << (entryToken + SEMI | testDict | testList | testList2 )
        #tempD = Dict(delimitedList(result,delim=SEMI)) \
        #       .ignore(cStyleComment).ignore('//' + restOfLine) \
        #       .parseString(expressions)

        #print tempD.dump()
        #f.close()

        return caseDict
        ##----------------------------------------------------------------------
