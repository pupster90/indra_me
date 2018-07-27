



# coding: utf-8

# <p></p><font size="7" face="courier" color="magenta">cx_assembler</font>

# This `.ipynb` file creates the class <font color=magenta>CxAssembler</font> that turns indra statements into a `.cx` file. Eventually this code is meant to replace indra's [cx_assmbler.py](
# https://github.com/sorgerlab/indra/blob/master/indra/assemblers/cx_assembler.py) file. Use the jupyter's *"export"* command to turn this notebook into a `.py` file. The [cx_assembler-DOCUMENTATION](http://indra.akintunde.org/notebooks/source/cx_assembler-DOCUMENTATION.ipynb) file gives a more thorough explanation about what this notebook does.

# ## <font color=red>Issues and Future Work :</font>

# * <font color=orange>Have statementCreator file call indraView to make code more concise</font>
# * <font color=red>If there is no name id's in _get_agent_alias we simply print "NA"</font>
# * <font color=red>Too many of the nodes are of type "other" check out the function,</font><font color=blue>_get_agent_type</font>
# * <font color=red>Figure out how to deal with the following warnings:</font>
# <img src="../data/Images/bel_warnings.png" width=700>
# 

# # Setup Notebook

# ## Setup Notebook Display

# In[1]:


from IPython.core.display import HTML, display


# In[3]:


# to reverse this feature restart kernel with "restart and clear output"
"""
display(HTML("<style> *{margin:0; padding:0;} html, body, \
             .container{margin:0;!important padding:0;!important} \
             .container { width:100% !important;}</style>"))
"""


# ## Load Stuff
# 
# We import the exact packages that are used by default from [cx_assmbler.py](
# http://35.202.250.158/edit/indra/indra/assemblers/cx_assembler.py)

# We add the path to the indra folder so that the commands of the following form work :
# ``` python
# from indra... 
# ```


import sys
sys.path.append("/root/Documents/indra")


# Import pacakges. These are the packages that the default 
# 

# In[5]:



from builtins import dict, str
import io
import re
import json
import logging
import itertools
from collections import OrderedDict
from indra.statements import *
from indra.databases import context_client, ndex_client, get_identifiers_url


# In[1]:


import warnings


# In[4]:


import ndex2


# Python 2
try:
    basestring
# Python 3
except:
    basestring = str

logger = logging.getLogger('cx_assembler')


# ndex2 (development edition) is used to create CX file

# In[5]:


# Code used for debugging stuff
#import inspect
#inspect.getfile(ndex2)





# # <font color="magenta">CxAssembler</font> Class
# This initializes the <font color="magenta">CxAssembler</font> class. Then in the later sections we add functions to this class

# ## Class Initialisation

# In[1]:


class CxAssembler(object):
    """This class assembles a CX network from a set of INDRA Statements.

    The CX format is an aspect oriented data mode for networks.
    The format is defined at http://www.home.ndexbio.org/data-model/.
    The CX format is the standard for NDEx and is compatible with
    CytoScape via the CyNDEx plugin.

    Parameters
    ----------
    stmts : Optional[list[indra.statements.Statement]]
        A list of INDRA Statements to be assembled.

    Attributes
    ----------
    statements : list[indra.statements.Statement]
        A list of INDRA Statements to be assembled.
    # network_name : str  <-- I moved this to the make_model function for consitency
        The name of the network to be assembled.
    cx : dict
        The structure of the CX network that is assembled.
    """
   
    def __init__( self, stmts ):
                
        # We create the CX file using the NDEX2 niceCX package
        self.myCX  = ndex2.NiceCXNetwork()
        
        # We save the save the statements setup some other helper variables
        self.statements = stmts
        self._existing_edges = {}
        


# In[60]:


print("Imported: CxAssembler")


# # <font color="blue">_add</font>  Functions
# 
# As the name implies, all functions that start with <font color="blue">_add</font> are used to add things to already prexisting objects. The <font color=gray>Add Statements</font> subsection consists of the functions that append additional statemens to the dataset. The <font color=red>Add CX Data</font> functions are used to turn the information from indra statements into a CX file. 
# 
# 
# 
# *<font color=gray>Future Work: Creating sub-objects within the <font color=magenta>CxAssembler</font> class  and divvying up the functions into the sub-objects would make this code more clear*
# 
# 

# ## <font color="gray"> Dredge Statements </font> 

# ### _add_modification 

# In[8]:


def _add_modification(self, stmt):
    if stmt.enz is None:
        return
    enz_id = self._add_node(stmt.enz)
    sub_id = self._add_node(stmt.sub)
    stmt_type = stmt.__class__.__name__
    self._add_edge(enz_id, sub_id, stmt_type, stmt)


# In[9]:


CxAssembler._add_modification = _add_modification


# ### _add_self_modification

# In[10]:


def _add_self_modification(self, stmt):
    enz_id = self._add_node(stmt.enz)
    stmt_type = stmt.__class__.__name__
    self._add_edge(enz_id, enz_id, stmt_type, stmt)


# In[11]:


CxAssembler._add_self_modification = _add_self_modification


# ### _add_complex

# In[12]:


def _add_complex(self, stmt):
    for m1, m2 in itertools.combinations(stmt.members, 2):
        m1_id = self._add_node(m1)
        m2_id = self._add_node(m2)
        self._add_edge(m1_id, m2_id, 'Complex', stmt)


# In[13]:


CxAssembler._add_complex = _add_complex


# ### _add_regulation

# In[14]:


def _add_regulation(self, stmt):
    if stmt.subj is None:
        return
    subj_id = self._add_node(stmt.subj)
    obj_id = self._add_node(stmt.obj)
    stmt_type = stmt.__class__.__name__
    self._add_edge(subj_id, obj_id, stmt_type, stmt)


# In[15]:


CxAssembler._add_regulation = _add_regulation


# ### _add_gef

# In[16]:


def _add_gef(self, stmt):
    gef_id = self._add_node(stmt.gef)
    ras_id = self._add_node(stmt.ras)
    stmt_type = stmt.__class__.__name__
    self._add_edge(gef_id, ras_id, stmt_type, stmt)


# In[17]:


CxAssembler._add_gef = _add_gef


# ### _add_gap

# In[18]:


def _add_gap(self, stmt):
    gap_id = self._add_node(stmt.gap)
    ras_id = self._add_node(stmt.ras)
    stmt_type = stmt.__class__.__name__
    self._add_edge(gap_id, ras_id, stmt_type, stmt)


# In[19]:


CxAssembler._add_gap = _add_gap


# ## <font color="red"> CX Aspects </font>

# ### _add_node

# This functions adds a node to the <font color=magenta>nodes</font> and also adds the node's extra parameters to <font color=magenta>node attributes</font>

# In[3]:


def _add_node(self, agent):
    
    # If node already there, just return the node name
    node = self.myCX.get_node( agent.name )
    if type(node) != type(None) :
        return agent.name
    
    # Add node
    alias = _get_agent_alias( agent )
    node_id = alias.pop()
    node_name = agent.name  # node_id.split(':')[-1]
    node = self.myCX.create_node(id=node_id, node_name= node_name , node_represents= node_id)
    
    # Add node attributes
    self.myCX.add_node_attribute( property_of= node, name='type', values= _get_agent_type(agent)  ) #<-- adds type
    if len(alias) != 0 :
        self.myCX.add_node_attribute( property_of= node, name='alias', values= alias, type='list_of_string'  )  #<-- adds alias
    return node


# In[16]:


CxAssembler._add_node = _add_node


# ### _add_edge

# In[10]:


def _add_edge(self, source, target, interaction, stmt):
    
    # if edge is already in list of existing edges, just return edge id
    edge_key = (source, target, interaction)
    try:
        edge_id = self._existing_edges[edge_key]
        return edge_id
    except KeyError:
        pass
    
    # add edge
    interaction = 'associates-with' if stmt.__class__.__name__ == 'Complex' else 'regulates'
    edge_id = self.myCX.create_edge( edge_source=source, edge_target=target, edge_interaction=interaction )
    self._existing_edges[edge_key] = edge_id  #<-- add edge id to list of existing edges
    self._add_edge_metadata(edge_id, stmt)    #<-- add edge attributes
    return edge_id


# In[67]:


CxAssembler._add_edge = _add_edge


# ### _add_edge_metadata

# In[9]:


def _add_edge_metadata(self, edge_id, stmt):
    
    # Add mechanism
    mechanism = _get_stmt_mechanism(stmt)
    if mechanism != "NA" :
        self.myCX.add_edge_attribute( property_of= edge_id, name='mechanism', values= mechanism )
    
    # Add polarity and contact
    self.myCX.add_edge_attribute( property_of= edge_id, name='polarity', values= _get_stmt_type(stmt)[1] ) #<-- polarity
    self.myCX.add_edge_attribute( property_of= edge_id, name='contact', values= stmt.evidence[0].epistemics["direct"] )  #<-- contact
    
    # Add Interaction
    #interaction = 'associates-with' if stmt.__class__.__name__ == 'Complex' else 'regulates'
    #self.myCX.add_edge_attribute( property_of= edge_id, name='interaction', values= interaction )  
    
    # Add the Citations for the Edge
    citations = _get_stmt_citations(stmt)
    if citations != "[]" :
        self.myCX.add_edge_attribute( property_of= edge_id, name='citations', values= citations, type='list_of_string' )
    
    # Add Evidence
    self.myCX.add_edge_attribute( property_of= edge_id, name='evidence', values=[stmt.evidence[0].text], type='list_of_string' )
    
    # For debugging, add indra statement, indra trext, and the indra JSON
    if self.debug ==  True :
        self.myCX.add_edge_attribute( property_of= edge_id, name='indra', values= stmt.__str__()  )
        self.myCX.add_edge_attribute( property_of= edge_id, name='indra_JSON', values= str(stmt.to_json())  )
        self.myCX.add_edge_attribute( property_of= edge_id, name='indra_number', values= stmt.number  )


# In[69]:


CxAssembler._add_edge_metadata = _add_edge_metadata


# # <font color="blue">Main</font> Functions

# ##  <font color="red">save_model</font>

# Save the assembled CX network in a file.

# In[34]:


def save_model(self, file_name='model.cx'):
    
    with open(file_name, 'wt') as fh:
        fh.write( self.cx_json )
        


# In[35]:


CxAssembler.save_model = save_model


# ## <font color="red">make_model</font>
# 

# <font color="red">make_model</font> adds the indra statements according which of the <font color="blue">_add</font> commands were specified. It is also where the user decides the <font color=green>network attribute</font> properties. <font color="red">make_model</font> returns the assembled CX network as a json string.
# 

# In[17]:


def make_model(self, name='indra_assembled', description='An Indra Auto-Curated network', version='1.0', debug=True ):
    
    """Assemble the CX network from the collected INDRA Statements.

    This method assembles a CX network from the set of INDRA Statements.
    The assembled network is set as the assembler's cx argument.

    Parameters
    ----------
    add_indra_json : Optional[bool]
        If True, the INDRA Statement JSON annotation is added to each
        edge in the network. Default: True

    Returns
    -------
    cx_str : str
        The json serialized CX model.
    """
    # save debug status, this is used later in hte "add" functions
    self.debug = debug
    
    # add @context to the CX file
    self.myCX.set_namespaces([{ "pubmed": "http://identifiers.org/pubmed/",
        "HGNC": "http://identifiers.org/hgnc/",
        "MGI": "http://www.informatics.jax.org/searchtool/Search.do?query=", 
        "RGD": "https://rgd.mcw.edu/rgdweb/elasticResults.html?category=Gene&species=Rat&cat1=General&sp1=&postCount=1&term=",
        "chebi": "http://identifiers.org/CHEBI:",
        "uniprot":"http://www.uniprot.org/uniprot/?query=",
        # Still have to figure out: (below are old db's not sure if actually used)
        "cas": "http://identifiers.org/cas/",
        "hprd": "http://identifiers.org/hprd/", 
        "KEGG Compound": "http://identifiers.org/kegg.compound/" }])

    # Add extra indra statements depending on settings
    for stmt in self.statements:
        if isinstance(stmt, Modification):
            self._add_modification(stmt)
        if isinstance(stmt, SelfModification):
            self._add_self_modification(stmt)
        elif isinstance(stmt, RegulateActivity) or isinstance(stmt, RegulateAmount):
            self._add_regulation(stmt)
        elif isinstance(stmt, Complex):
            self._add_complex(stmt)
        elif isinstance(stmt, Gef):
            self._add_gef(stmt)
        elif isinstance(stmt, Gap):
            self._add_gap(stmt)
            
    # Add network Attributes
    self.myCX.set_name( name )
    self.myCX.add_network_attribute(name='description', values= description  )#, type=ATTRIBUTE_DATA_TYPE.STRING) #<- removed type variable
    self.myCX.add_network_attribute(name='version', values= version )#, type=ATTRIBUTE_DATA_TYPE.STRING)  #<- removed type variable
    
    self.cx_json = json.dumps( self.myCX.to_cx() ,indent=2)
    return self.cx_json


# In[18]:


CxAssembler.make_model = make_model


# ## <font color="red">upload_model</font>

# Creates a new NDEx network of the assembled CX model.
# 
#     To upload the assembled CX model to NDEx, you need to have
#     a registered account on NDEx (http://ndexbio.org/) and have
#     the `ndex` python package installed. The uploaded network
#     is private by default.

# In[79]:


def upload_model(self, account, password,  url="http://www.ndexbio.org"):
    
    try:
        my_ndex= ndex2.client.Ndex2( url, account, password)
        my_ndex.update_status()
        print("Success.  Please continue.")
    except Exception as inst:
        print("Could not access account %s with password %s" % (account, password))
        print(inst.args)
    upload_message  = self.myCX.upload_to( url, account, password)
    print(upload_message)


# In[41]:


CxAssembler.upload_model = upload_model


# # <font color="blue">_get</font> functions 

# ## <font color=magenta>Node Attributes</font>

# ### _get_agent_type

# In[2]:


def _get_agent_type(agent):

    
    mgi_id = agent.db_refs.get('MGI')
    rgd_id = agent.db_refs.get('RGD')
    hgnc_id = agent.db_refs.get('HGNC')
    if mgi_id or rgd_id or hgnc_id:
        return 'gene'
        
    uniprot_id = agent.db_refs.get('UP')
    uniprot2_id = agent.db_refs.get('EGID')
    if uniprot_id or uniprot2_id:
        return 'protein'
        
    pfam_id = agent.db_refs.get('PF')
    fa_id = agent.db_refs.get('FA')
    be_id = agent.db_refs.get('FPLX')
    if pfam_id or fa_id or be_id:
        return 'protein_family'
    
    chebi_id = agent.db_refs.get('CHEBI')
    pubchem_id = agent.db_refs.get('PUBCHEM')
    schem_id = agent.db_refs.get('SCHEM')
    if chebi_id or pubchem_id or schem_id :
        return 'chemical'

    #go_id = agent.db_refs.get('GO')   
    #if go_id:
    #    return 'biological process'
    
    # Some other type of id
    return 'other'


# In[49]:


CxAssembler._get_agent_type = _get_agent_type


# ### _get_agent_alias

# <font color=red>Warning: in current code if there is no name id's we simply print "NA"</font>

# In[11]:


# helper function

# This code add's the alias'es for a node (if they exist)
def _get_agent_alias(agent):
    alias = []
    #print(" ");print(agent.db_refs)
    
    for db_name, db_ids in agent.db_refs.items():
        if not db_ids:
            logger.warning('Missing db_id for %s' % agent)
            continue
        elif isinstance(db_ids, int):
            db_id = str(db_ids)
        elif isinstance(db_ids, basestring):
            db_id = db_ids.split(':')[-1] #<-- get only id if there's an namespace:id format 
        else:
            db_id = db_ids[0]
        
        db_name_map = {
            'UP': 'uniprot knowledgebase', 'EGID':'uniprot knowledgebase',
            'PUBCHEM': 'PubChem', 'HGNC':'HGNC' ,
            'IP': 'InterPro', 'NXPFA': 'NextProtFamily',
            'PF': 'Pfam', 'CHEBI': 'chebi', 
            'MGI': 'MGI', 'RGD': 'RGD' }
        
        name = db_name_map.get(db_name)
        if not name:
            warnings.warn("NAME: "+db_name+" , not found" )
            name = db_name
        alias.append( name+":"+db_id ) 
    #print(alias)
    
    if len(alias)==0:
        warnings.warn("Alias length of 0")
        return ["NA:"+agent.name]  #raise ValueError('indra statement contains no ID!')
    return alias


# In[51]:


CxAssembler._get_agent_alias = _get_agent_alias


# ## <font color=magenta>Edge Attributes</font>

# ### _get_stmt_type
# 
# This could be a useful spot for  creating the specified `type` attributes in the [CX Documenation](https://docs.google.com/document/d/13ZKcFBH-E5oiJP2D5zrdFqxLlS9yGtGiiX5Lj92g4EU/edit#heading=h.1t5lf4irpgyj)

# In[46]:


def _get_stmt_type(stmt):
    if isinstance(stmt, AddModification):
        edge_type = 'Modification'
        edge_polarity = 1
    elif isinstance(stmt, RemoveModification):
        edge_type = 'Modification'
        edge_polarity = -1
    elif isinstance(stmt, SelfModification):
        edge_type = 'SelfModification'
        edge_polarity = 1
    elif isinstance(stmt, Complex):
        edge_type = 'Complex'
        edge_polarity = 0
    elif isinstance(stmt, Activation):
        edge_type = 'Activation'
        edge_polarity = 1
    elif isinstance(stmt, Inhibition):
        edge_type = 'Inhibition'
        edge_polarity = -1
    elif isinstance(stmt, DecreaseAmount):
        edge_type = 'DecreaseAmount'
        edge_polarity = -1
    elif isinstance(stmt, IncreaseAmount):
        edge_type = 'IncreaseAmount'
        edge_polarity = 1
    elif isinstance(stmt, Gef):
        edge_type = 'Gef'
        edge_polarity = 1
    elif isinstance(stmt, Gap):
        edge_type = 'Gap'
        edge_polarity = -1
    else:
        edge_type = stmt.__class__.__str__()
        edge_polarity = 'none'
    return edge_type, edge_polarity


# In[47]:


CxAssembler._get_stmt_type = _get_stmt_type


# ### _get_agent_citations

# In[52]:



def _get_stmt_citations(stmt):
    pmids = [e.pmid for e in stmt.evidence if e.pmid]
    pmids_added = []
    for pmid in pmids:
        pmid_txt = None
        if re.match('[0-9]+', pmid):
            pmid_txt = 'pubmed:' + pmid
            if pmid_txt not in pmids_added:
                pmids_added.append(pmid_txt)
    return pmids_added


# In[53]:


CxAssembler._get_stmt_citations = _get_stmt_citations


# ### _get_stmt_mechanism
# 
# maps the indra mechanisms into the desired CX mechanism

# In[14]:


def _get_stmt_mechanism( stmt ):
    pre_map_name = stmt.__class__.__name__ if stmt.__class__.__name__ != "Activation" else stmt.obj_activity
    
    if pre_map_name in ["Acetylation","Farnesylation","Ubiquitination","Phosphorylation","Methylation"] :
        return pre_map_name.lower()
    
    if pre_map_name in ["Gap","Gef","GtpActivation","gtpbound"]:
        return "gtp bound"
    
    if pre_map_name == "transcription":
        return "transcriptional regulation"
    
    if pre_map_name == "phosphatase" :
        return "dephosphorylation"
    
    if pre_map_name == "kinase" :
        return "phosphorylation"
    
    if pre_map_name == "Complex":
        return "binding"
    
    if pre_map_name == "catalytic":
        return "small molecule catalysis reaction"
    
    return "NA"
    


# In[12]:


CxAssembler._get_stmt_mechanism = _get_stmt_mechanism

