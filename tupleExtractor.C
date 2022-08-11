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
  TFile *f = new TFile("LQD_tuple.root", "RECREATE");

  // Set Input Filenames in the string array
  std::vector<std::string> inputFileNames;
  inputFileNames.push_back("DAOD_TRUTH0.test.pool.truth0.root");

  // Initialize xAOD stuff
  auto xaodEvent = new xAOD::TEvent(xAOD::TEvent::kClassAccess);
  TFile *inFile = 0;

  // Provide Tree name for the ntuples/branches to be written out
  TTree *tree = new TTree("LQDTruthTuple", "Truth0 Information");

  // Declare data types; must match the data type being pulled out of TRUTH0 DAOD
  int eventNumber;
  std::vector<float> DV_R;
  std::vector<float> n1_lifetime;

  std::vector<int> nn1;
  std::vector<int> ng;

  std::vector<int> njets;
  std::vector<int> nbjets;

  std::vector<int> nleptons;
  std::vector<int> nelectrons;
  std::vector<int> nmuons;

  std::vector<float> n1pt;
  std::vector<float> gpt;

  std::vector<float> jpt;
  std::vector<float> bjpt;

  std::vector<float> leppt;
  std::vector<float> ept;
  std::vector<float> mupt;

  // Initialize the branches
  tree->Branch("EventNumber", &eventNumber);
  tree->Branch("DV_R", &DV_R);
  tree->Branch("n1_lifetime", &n1_lifetime);

  tree->Branch("Nn1", &nn1);
  tree->Branch("Ng", &ng);
  tree->Branch("n1pt", &n1pt);
  tree->Branch("gpt", &gpt);

  tree->Branch("Njets", &njets);
  tree->Branch("Nbjets", &nbjets);

  tree->Branch("Nleptons", &nleptons);
  tree->Branch("Nelectrons", &nelectrons);
  tree->Branch("Nmuons", &nmuons);

  tree->Branch("jpt", &jpt);
  tree->Branch("bjpt", &bjpt);

  tree->Branch("leppt", &leppt);
  tree->Branch("ept", &ept);
  tree->Branch("mupt", &mupt);

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

      int nugluinos = 0;
      int nuneutralinos = 0;

      int nubjets = 0;

      int nuleptons = 0;
      int nuelectrons = 0;
      int numuons = 0;

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
          nugluinos++;
          gpt.push_back(1e-3 * SP->pt());
        }

        if (fabs(SP->pdgId()) == 5 && SP->status() == 1)
        {
          nubjets++;
          bjpt.push_back(1e-3 * SP->pt());
          if (!DEBUG)
            std::cout << "The b-jet Pt is:" << SP->pt() << std::endl;
        }

        if (fabs(SP->pdgId()) >= 11 && fabs(SP->pdgId()) <= 18 && SP->status() == 1)
        {
          nuleptons++;
          leppt.push_back(1e-3 * SP->pt());
        }
        if (fabs(SP->pdgId()) == 11 && SP->status() == 1)
        {
          nuelectrons++;
          ept.push_back(1e-3 * SP->pt());
          if (!DEBUG)
            std::cout << "The Electron Pt is:" << SP->pt() << std::endl;
        }
        if (fabs(SP->pdgId()) == 13 && SP->status() == 1)
        {
          numuons++;
          mupt.push_back(1e-3 * SP->pt());
        }

        if (fabs(SP->pdgId()) == 1000022 && SP->status() == 22)
        {
          hard_int.push_back(tp);
          nuneutralinos++;
          n1pt.push_back(1e-3 * SP->pt());

          if (DEBUG)
            std::cout << "pdgID: " << SP->pdgId() << ", mass: " << SP->m() / 1000. << ", decays? " << SP->hasDecayVtx() << ", status: " << SP->status() << std::endl;
          if (DEBUG)
            std::cout << "Production vertex: (" << SP->prodVtx()->x() << ", " << SP->prodVtx()->y() << ", " << SP->prodVtx()->z() << ")" << std::endl;
          if (DEBUG)
            std::cout << "Decay vertex: (" << SP->decayVtx()->x() << ", " << SP->decayVtx()->y() << ", " << SP->decayVtx()->z() << ")" << std::endl;
          float decayR = sqrt(pow(((SP->prodVtx()->x()) - (SP->decayVtx()->x())), 2.) + pow(((SP->prodVtx()->y()) - (SP->decayVtx()->y())), 2.) + pow(((SP->prodVtx()->z()) - (SP->decayVtx()->x())), 2.));
          DV_R.push_back(0.1 * decayR);

          // Is in seconds, because Physics Short was used
          float lifetimeLab = (SP->decayVtx()->v4().Vect()).Mag() / (SP->p4().Beta() * SP->p4().Gamma() * TMath::C()) / 1000.;
          if (DEBUG)
            std::cout << "lifetime: " << lifetimeLab << std::endl
                      << std::endl;
          if (lifetimeLab > 0.)
            n1_lifetime.push_back(lifetimeLab);
        }
      } // end of loop over truth particles

      // Fill number of respective particles in each event
      nn1.push_back(nuneutralinos);
      ng.push_back(nugluinos);
      nbjets.push_back(nubjets);
      nleptons.push_back(nuleptons);
      nelectrons.push_back(nuelectrons);
      nmuons.push_back(numuons);

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
      n1_lifetime.clear();
    } // end of entries loop
    inFile->Close();
  } // end of filenames loop
  f->cd();
  f->WriteObject(tree, "LQDTruthTuple");
  tree->SetDirectory(f);
  tree->Write();
}
