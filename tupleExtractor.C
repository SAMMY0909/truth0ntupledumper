/**
 * @brief This code extracts ntuples out of TRUTH0 DAODs
 * @authors Emily Anne Thompson, Soumyananda Goswami
 * Initialize environment before running, on any machine that has cvmfs access:
 * `setupATLAS -q`
 * `asetup AnalysisBase,21.2.125`
 * To run:
 * `root -l -b -q tupleExtractor.C`
 */
bool wasProducedFromSMquarks(const xAOD::TruthParticle *truthParticle)
{
  bool fromSMQuarks = true;
  if (!truthParticle->hasProdVtx())
    return false;
  for (size_t i = 0; i < truthParticle->prodVtx()->nIncomingParticles(); i++)
  {
    if (fabs(truthParticle->prodVtx()->incomingParticle(i)->pdgId()) > 9)
      fromSMQuarks = false;
  }
  return fromSMQuarks;
}

void fillChildMap(std::map<std::pair<int, int>, const xAOD::TruthParticle *> *childMap, const xAOD::TruthParticle *parent)
{

  auto parent_pdgId = parent->pdgId();
  if (!parent->decayVtx())
  {
  } // do nothing
  else
  {
    for (size_t i = 0; i < parent->decayVtx()->nOutgoingParticles(); i++)
    {

      const xAOD::TruthParticle *child = parent->decayVtx()->outgoingParticle(i);

      if (child)
      {
        // if (child->status()==1 ){
        // if (fabs(child->pdgId())==1000022 || !(fabs(child->pdgId())>1e6 && fabs(child->pdgId())<3e6)){
        if (fabs(child->pdgId()) != parent_pdgId)
        {
          childMap->insert(std::pair<std::pair<int, int>, const xAOD::TruthParticle *>(std::make_pair(child->pdgId(), child->barcode()), child));
        }
        else
        {
          fillChildMap(childMap, child);
        }
      }
    }
  }
  return;
}

void tupleExtractor()
{
  // Set Debug state
  bool DEBUG = false;

  // Provide Output File Name
  TFile *f = new TFile("LQD_tuple.root", "NEW");

  // Set Input Filenames in the string array
  std::vector<std::string> inputFileNames;
  inputFileNames.push_back("DAOD_TRUTH0.test.pool.truth0.root");

  // Initialize xAOD stuff
  auto xaodEvent = new xAOD::TEvent(xAOD::TEvent::kClassAccess);
  TFile *inFile = 0;

  // Provide Tree name for the ntuples/branches to be written out
  TTree *tree = new TTree("LQDTruthTuple", "Truth0 Information");

  // Declare data types; must match the data type being pulled out of TRUTH0 DAOD
  int m_eventNumber;
  std::vector<float> m_DV_R;
  std::vector<float> m_n1_lifetime;

  std::vector<int> m_nn1;
  std::vector<int> m_ng;

  std::vector<int> m_njets;
  std::vector<int> m_nbjets;

  std::vector<int> m_nleptons;
  std::vector<int> m_nelectrons;
  std::vector<int> m_nmuons;

  std::vector<float> m_jpt;
  std::vector<float> m_bjpt;

  std::vector<float> m_leppt;
  std::vector<float> m_ept;
  std::vector<float> m_mupt;

  // Initialize the branches
  tree->Branch("EventNumber", &m_eventNumber);
  tree->Branch("DV_R", &m_DV_R);
  tree->Branch("n1_lifetime", &m_n1_lifetime);

  tree->Branch("Nn1", &m_nn1);
  tree->Branch("Ng", &m_ng);

  tree->Branch("Njets", &m_njets);
  tree->Branch("Nbjets", &m_nbjets);

  tree->Branch("Nleptons", &m_nleptons);
  tree->Branch("Nelectrons", &m_nelectrons);
  tree->Branch("Nmuons", &m_nmuons);

  tree->Branch("jpt", &m_jpt);
  tree->Branch("bjpt", &m_bjpt);

  tree->Branch("leppt", &m_leppt);
  tree->Branch("ept", &m_ept);
  tree->Branch("mupt", &m_mupt);

  // Loop over filenames
  for (const auto &inFileName : inputFileNames)
  {
    // Delete any previous input files initialized in the inFile array
    delete inFile;

    // Open DAOD TRUTH0 file as READ only
    inFile = TFile::Open(inFileName.c_str(), "READ");
    if (!xaodEvent->readFrom(inFile).isSuccess())
    {
      throw std::runtime_error("Could not connect TEvent to file!");
    }

    // Get Number of Events, must equal what you specified in the MC production, before cuts
    Long64_t numEntries = xaodEvent->getEntries();
    std::cout << "Num Event Entries=" << numEntries << std::endl;

    // if (DEBUG) numEntries = 20; //This is if you need fewer events to catch errors.
    // Loop over events
    for (Long64_t index = 0; index < numEntries; index++)
    {

      // Get n-th event
      Long64_t entry = xaodEvent->getEntry(index);
      if (entry < 0)
      {
        std::cout << "Entry less than 0!" << std::endl;
      }
      if (DEBUG)
        std::cout << "================= New event =====================" << std::endl
                  << std::endl;

      // Get basic event info
      const xAOD::EventInfo *eventInfo = 0;
      if (!xaodEvent->retrieve(eventInfo, "EventInfo").isSuccess())
      {
        throw std::runtime_error("Cannot read Event Info");
      }

      // Get truth particles
      /**
       * For Truth Particle Container, Please see: https://ucatlas.github.io/RootCoreDocumentation/2.4.28/dd/dc2/classxAOD_1_1TruthParticle__v1.html
       * For Truth Vertex Container, Please see: https://ucatlas.github.io/RootCoreDocumentation/2.4.28/d8/dfa/classxAOD_1_1TruthVertex__v1.html
       */
      const xAOD::TruthParticleContainer *truthparticles = 0;
      if (xaodEvent->contains<xAOD::TruthParticleContainer>("TruthParticles"))
      {
        if (!xaodEvent->retrieve(truthparticles, "TruthParticles").isSuccess())
        {
          throw std::runtime_error("Could not retrieve truth particles");
        }
      }
      // if (DEBUG)
      std::cout << "Number of truth particles in this event are: " << truthparticles->size() << std::endl;

      vector<size_t> hard_int = {};

      int num_gluinos = 0;
      int num_neutralinos = 0;

      int num_bjets = 0;

      int num_leptons = 0;
      int num_electrons = 0;
      int num_muons = 0;

      // Check the particles produced directly from hard interactions
      if (DEBUG)
        std::cout << "Particles produced in hard interaction: " << std::endl;
      for (size_t tp = 0; tp < truthparticles->size(); tp++)
      {
        const auto SP = truthparticles->at(tp); // Get Selected truth particle at the given iterator index

        /**
         * The following lines of code select particles by their PDGID
         * Please see: https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
         * The particles are counted by their respective type per event
         * Their pts are also filled up, where required
         */
        if ((fabs(SP->pdgId()) == 1000021) && SP->status() == 22)
        {
          num_gluinos++;
        }

        if ((fabs(SP->pdgId()) == 5) && SP->status() == 22)
        {
          num_bjets++;
          m_bjpt.push_back(SP->pt());
        }

        if (((fabs(SP->pdgId()) >= 11) && SP->status() == 22) || ((fabs(SP->pdgId()) <= 18) && SP->status() == 22))
        {
          num_leptons++;
          m_leppt.push_back(SP->pt());
        }
        if ((fabs(SP->pdgId()) == 11) && SP->status() == 22)
        {
          num_electrons++;
          m_ept.push_back(SP->pt());
        }
        if ((fabs(SP->pdgId()) == 13) && SP->status() == 22)
        {
          num_muons++;
          m_mupt.push_back(SP->pt());
        }

        if ((fabs(SP->pdgId()) == 1000022) && SP->status() == 22)
        {
          hard_int.push_back(tp);
          num_neutralinos++;

          if (DEBUG)
            std::cout << "pdgID: " << SP->pdgId() << ", mass: " << SP->m() / 1000. << ", decays? " << SP->hasDecayVtx() << ", status: " << SP->status() << std::endl;
          if (DEBUG)
            std::cout << "Production vertex: (" << SP->prodVtx()->x() << ", " << SP->prodVtx()->y() << ", " << SP->prodVtx()->z() << ")" << std::endl;
          if (DEBUG)
            std::cout << "Decay vertex: (" << SP->decayVtx()->x() << ", " << SP->decayVtx()->y() << ", " << SP->decayVtx()->z() << ")" << std::endl;
          float decayR = sqrt(pow(((SP->prodVtx()->x()) - (SP->decayVtx()->x())), 2.) + pow(((SP->prodVtx()->y()) - (SP->decayVtx()->y())), 2.) + pow(((SP->prodVtx()->z()) - (SP->decayVtx()->x())), 2.));
          m_DV_R.push_back(decayR);

          // Calculated in E-3 s. Multiply by 1000 to get to seconds
          float lifetimeLab = (SP->decayVtx()->v4().Vect()).Mag() / (SP->p4().Beta() * SP->p4().Gamma() * TMath::C()) / 1000.;
          if (DEBUG)
            std::cout << "lifetime: " << lifetimeLab << std::endl
                      << std::endl;
          if (lifetimeLab > 0.)
            m_n1_lifetime.push_back(lifetimeLab);
        }
      } // end of loop over truth particles

      // Fill number of respective particles in each event
      m_nn1.push_back(num_neutralinos);
      m_ng.push_back(num_gluinos);
      m_nbjets.push_back(num_bjets);
      m_nleptons.push_back(num_leptons);
      m_nelectrons.push_back(num_electrons);
      m_nmuons.push_back(num_muons);

      if (DEBUG)
        std::cout << std::endl
                  << "Stable particles from hard interaction decays: " << std::endl;

      // The : stands for fancy C++ iterator and loops (in this case) over all hard interaction particles
      for (auto tp : hard_int)
      {
        const auto SP = truthparticles->at(tp);

        std::map<std::pair<int, int>, const xAOD::TruthParticle *> childMap;
        fillChildMap(&childMap, SP);
        std::vector<const xAOD::TruthParticle *> uniqueChildren;
        typedef std::map<std::pair<int, int>, const xAOD::TruthParticle_v1 *>::const_iterator MapIterator;
        for (MapIterator iter = childMap.begin(); iter != childMap.end(); iter++)
        {
          uniqueChildren.push_back(iter->second);
          if (DEBUG)
            std::cout << "Daughter pdgId: " << iter->second->pdgId() << ", mass: " << iter->second->m() / 1000. << ", pt: " << iter->second->pt() / 1000. << std::endl;
        }
      }
      tree->Fill();
      m_n1_lifetime.clear();
    } // end of entries loop
    inFile->Close();
  } // end of filenames loop
  f->cd();
  f->WriteObject(tree, "LQDTruthTuple");
  tree->SetDirectory(f);
  tree->Write();
}
