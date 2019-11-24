#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"

#include <vector>
#include <algorithm>

using namespace std;
using namespace edm;

namespace flashgg {

    class SinglePhotonViewProducer : public EDProducer
    {

    public:
        SinglePhotonViewProducer( const ParameterSet & );
    private:
        SinglePhotonView MakeRecHitMap(SinglePhotonView, Handle<EcalRecHitCollection>, Handle<EcalRecHitCollection>, const CaloSubdetectorGeometry*, const CaloSubdetectorGeometry*);
        void produce( Event &, const EventSetup & ) override;

        EDGetTokenT<View<DiPhotonCandidate> > diPhotonToken_;
        int maxCandidates_;
        edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > EB_reducedEcalRecHitsToken_;
        edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > EE_reducedEcalRecHitsToken_;

    };

    SinglePhotonViewProducer::SinglePhotonViewProducer( const ParameterSet &iConfig ) :
        diPhotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
        maxCandidates_( iConfig.getParameter<int>( "maxCandidates" ) ),
        EB_reducedEcalRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("EBreducedEcalRecHits"))),
        EE_reducedEcalRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("EEreducedEcalRecHits")))
    {
        produces<vector<SinglePhotonView> >();
    }

    SinglePhotonView SinglePhotonViewProducer::MakeRecHitMap( SinglePhotonView spv, Handle<EcalRecHitCollection> EBReducedRecHits, Handle<EcalRecHitCollection> EEReducedRecHits,  const CaloSubdetectorGeometry* _ebGeom, const CaloSubdetectorGeometry* _eeGeom)
    {
        std::vector<float> DOF1s;
        std::vector<float> DOF2s;
        std::vector<float> DOF3s;
        std::vector<float> recHits; 
        
        unsigned int maxnHits = 100;
        float dummyValue = -9999; 
        for (unsigned int i = 0; i < maxnHits; ++i) {
            DOF1s.push_back(dummyValue);
            DOF2s.push_back(dummyValue);
            DOF3s.push_back(dummyValue);
            recHits.push_back(dummyValue);
        }

        int pho_DOF1; // ieta (ix) for EB (EE) 
        int pho_DOF2; // iphi (iy) for EB (EE) 
        // int pho_DOF3; // iz. 0 (+/-1) for EB (EE) 
        int hitNum = 0;

        double rechitE;
        int rechitDOF1;
        int rechitDOF2;
        int rechitDOF3;

        // want to save all rechits in 10x10 window (for now. later make bigger) around photon eta and phi 

        double pho_eta = spv.originalPhoton()->superCluster()->eta();
        // double pho_phi = spv.originalPhoton()->phi(); 

        float pho_x = spv.originalPhoton()->superCluster()->position().x();   
        float pho_y = spv.originalPhoton()->superCluster()->position().y(); 
        float pho_z = spv.originalPhoton()->superCluster()->position().z(); 

        // cout << "pho eta = " << pho_eta << endl;
        // cout << "pho phi = " << pho_phi << endl;

        // cout << "pho x = " << pho_x << endl;
        // cout << "pho y = " << pho_y << endl;
        // cout << "pho z = " << pho_z << endl;

        // EB photon 
        if (abs(pho_eta) < 1.479){

            EBDetId pho_eb_id(_ebGeom->getClosestCell(GlobalPoint(pho_x,pho_y,pho_z)));
            pho_DOF1 = pho_eb_id.ieta();
            pho_DOF2 = pho_eb_id.iphi();
            // pho_DOF3 = 0;

            // cout << "pho ieta = " << pho_ieta << endl;
            // cout << "pho iphi = " << pho_iphi << endl;

            hitNum = 0;

            for (const auto & obj : *EBReducedRecHits) {

                rechitE = obj.energy();
                const EBDetId& EBhitId = obj.detid(); 
                rechitDOF1 = EBhitId.ieta();
                rechitDOF2 = EBhitId.iphi();
                rechitDOF3 = 0;

                // Only save the rec hit eta, phi, energy if its eta and phi are within 5 of the photon ieta iphi 

                if ( (abs((double)rechitDOF1 - (double)pho_DOF1) < 5) && (abs((double)rechitDOF2 - (double)pho_DOF2) < 5) ){
                    // cout << "rechitiEta = " << rechitiEta << endl;
                    // cout << "rechitiPhi = " << rechitiPhi << endl;
                    // cout << "rechitE = " << rechitE << endl;

                    DOF1s[hitNum] = (float)rechitDOF1;
                    DOF2s[hitNum] = (float)rechitDOF2;
                    DOF3s[hitNum] = (float)rechitDOF3;
                    recHits[hitNum] = (float)rechitE;

                    hitNum += 1;
                }


            }

        }

        // EE photon 
        else{

            EEDetId pho_ee_id(_eeGeom->getClosestCell(GlobalPoint(pho_x,pho_y,pho_z)));
            pho_DOF1 = pho_ee_id.ix();
            pho_DOF2 = pho_ee_id.iy();

            // if (pho_z < 0) pho_DOF3 = -1;
            // else if (pho_z > 0) pho_DOF3 = 1;
            // else{
            //     cout << "EE photon has z position = 0. This is probably a problem." << endl;
            // }

            // cout << "pho ieta = " << pho_ieta << endl;
            // cout << "pho iphi = " << pho_iphi << endl;

            hitNum = 0;

            for (const auto & obj : *EEReducedRecHits) {

                rechitE = obj.energy();
                const EEDetId& EEhitId = obj.detid(); 
                rechitDOF1 = EEhitId.ix();
                rechitDOF2 = EEhitId.iy(); 
                rechitDOF3 = EEhitId.zside();

                // Only save the rec hit eta, phi, energy if its eta and phi are within 5 of the photon ieta iphi 

                if ( (abs((double)rechitDOF1 - (double)pho_DOF1) < 5) && (abs((double)rechitDOF2 - (double)pho_DOF2) < 5) ){
                    // cout << "rechitiEta = " << rechitiEta << endl;
                    // cout << "rechitiPhi = " << rechitiPhi << endl;
                    // cout << "rechitE = " << rechitE << endl;

                    DOF1s[hitNum] = (float)rechitDOF1;
                    DOF2s[hitNum] = (float)rechitDOF2;
                    DOF3s[hitNum] = (float)rechitDOF3;
                    recHits[hitNum] = (float)rechitE;

                    hitNum += 1;
                }

            }

        }

        spv.SetDOF1s(DOF1s);
        spv.SetDOF2s(DOF2s);
        spv.SetDOF3s(DOF3s);
        spv.SetRecHits(recHits);

        DOF1s.clear();
        DOF2s.clear();
        DOF3s.clear();
        recHits.clear();

        return spv; 
    }    

    void SinglePhotonViewProducer::produce( Event &evt, const EventSetup &iSetup )
    {

        //calo geometry
        edm::ESHandle<CaloGeometry> caloGeometry;
        iSetup.get<CaloGeometryRecord>().get(caloGeometry);
        // const CaloGeometry *geometry = caloGeometry.product();
        // _ebGeom = caloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
        const CaloSubdetectorGeometry* _ebGeom; 
        const CaloSubdetectorGeometry* _eeGeom; 
        _ebGeom = caloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
        _eeGeom = caloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);

        Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
        evt.getByToken( diPhotonToken_, diPhotons );

        std::unique_ptr<vector<SinglePhotonView> > photonViews( new vector<SinglePhotonView> );

        Handle<EcalRecHitCollection>   EBReducedRecHits;
        Handle<EcalRecHitCollection>   EEReducedRecHits;
        // edm::Handle<edm::View<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > > EBReducedRecHits;
        evt.getByToken(EB_reducedEcalRecHitsToken_, EBReducedRecHits);
        evt.getByToken(EE_reducedEcalRecHitsToken_, EEReducedRecHits);


        int nCand = maxCandidates_;
        //for(auto & dipho : diPhotons) {
        for( unsigned int i = 0 ; i < diPhotons->size(); i++ ) {

            // For updated PhotonID study
            // Want to match leading and subleading photons with an array of rechits 
            // in order to have rec hit images for prompt and fake photons 

            SinglePhotonView lView = *diPhotons->ptrAt( i )->leadingView();
            SinglePhotonView slView = *diPhotons->ptrAt( i )->subLeadingView();
            
            lView = MakeRecHitMap(lView, EBReducedRecHits, EEReducedRecHits, _ebGeom, _eeGeom);
            slView = MakeRecHitMap(slView, EBReducedRecHits, EEReducedRecHits, _ebGeom, _eeGeom);          

            // Should push back after saving nearby rechit info 
            photonViews->push_back(lView);
            photonViews->push_back(slView);

            if( --nCand == 0 ) { break; }
        }

        //// if( photonViews->size() != 0 ) {
        //// 	cout << photonViews->size() << endl;
        //// }
        evt.put( std::move( photonViews ) );

    }
}

typedef flashgg::SinglePhotonViewProducer FlashggSinglePhotonViewProducer;
DEFINE_FWK_MODULE( FlashggSinglePhotonViewProducer );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

