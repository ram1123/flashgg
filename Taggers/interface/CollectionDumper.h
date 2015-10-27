#ifndef flashgg_CollectionDumper_h
#define flashgg_CollectionDumper_h

#include <map>
#include <string>

#include "TH1.h"
#include "TTree.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "PhysicsTools/UtilAlgos/interface/BasicAnalyzer.h"
/// #include "PhysicsTools/FWLite/interface/ScannerHelpers.h"

#include "CategoryDumper.h"

#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "flashgg/Taggers/interface/StringHelpers.h"

#include "RooWorkspace.h"
#include "RooMsgService.h"

#include "flashgg/MicroAOD/interface/CutBasedClassifier.h"
#include "flashgg/MicroAOD/interface/ClassNameClassifier.h"
#include "flashgg/MicroAOD/interface/CutAndClassBasedClassifier.h"
#include "flashgg/Taggers/interface/GlobalVariablesDumper.h"
#include "flashgg/DataFormats/interface/PDFWeightObject.h"

#include "FWCore/Utilities/interface/Exception.h"

/**
   \class CollectionDumper
   \brief Example class that can be used to analyze pat::Photons both within FWLite and within the full framework

   This is an example for keeping classes that can be used both within FWLite and within the full
   framework. The class is derived from the BasicAnalyzer base class, which is an interface for
   the two wrapper classes EDAnalyzerWrapper and FWLiteAnalyzerWrapper. You can fin more information
   on this on WorkBookFWLiteExamples#ExampleFive.
*/

namespace flashgg {


    template <class T>
    class TrivialClassifier
    {
    public:
        TrivialClassifier( const edm::ParameterSet &cfg ) {}

        std::pair<std::string, int> operator()( const T &obj ) const { return std::make_pair( "", 0 ); }
    };

    template<class CollectionT, class CandidateT = typename CollectionT::value_type, class ClassifierT = TrivialClassifier<CandidateT> >
    class CollectionDumper : public edm::BasicAnalyzer
    {

    public:
        typedef CollectionT collection_type;
        typedef CandidateT candidate_type;
        typedef StringObjectFunction<CandidateT, true> function_type;
        typedef CategoryDumper<function_type, candidate_type> dumper_type;
        typedef ClassifierT classifier_type;
        // typedef std::pair<std::string, std::string> KeyT;
        typedef std::string KeyT;

        /// default constructor
        CollectionDumper( const edm::ParameterSet &cfg, TFileDirectory &fs );
        CollectionDumper( const edm::ParameterSet &cfg, TFileDirectory &fs, const edm::ConsumesCollector &cc ) : CollectionDumper( cfg, fs ) {};
        /// default destructor
        /// ~CollectionDumper();
        /// everything that needs to be done before the event loop
        void beginJob();
        /// everything that needs to be done after the event loop
        void endJob();
        /// everything that needs to be done during the event loop
        void analyze( const edm::EventBase &event );


    protected:
        double eventWeight( const edm::EventBase &event );
        vector<double> pdfWeights( const edm::EventBase &event );

        /// float eventWeight(const edm::EventBase& event);
        edm::InputTag src_, genInfo_, pdfWeightToken_;
        std::string processId_;
        double lumiWeight_;
        int maxCandPerEvent_;
        double sqrtS_;
        double intLumi_;
        std::string nameTemplate_;

        bool dumpTrees_;
        bool dumpWorkspace_;
        std::string workspaceName_;
        bool dumpHistos_, dumpGlobalVariables_;

        classifier_type classifier_;
        std::map< KeyT, bool> hasSubcat_;
        bool throwOnUnclassified_;
        
        // event weight
        float weight_;
        vector<double> pdfWeights_;
        int pdfWeightSize_;
        bool pdfWeightHistosBooked_;
        bool dumpPdfWeights_;
        int nPdfWeights_;

        //std::map<std::string, std::vector<dumper_type> > dumpers_; FIXME template key
        std::map< KeyT, std::vector<dumper_type> > dumpers_;
        RooWorkspace *ws_;
        /// TTree * bookTree(const std::string & name, TFileDirectory& fs);
        /// void fillTreeBranches(const flashgg::Photon & pho)

        GlobalVariablesDumper *globalVarsDumper_;
    };

    template<class C, class T, class U>
    CollectionDumper<C, T, U>::CollectionDumper( const edm::ParameterSet &cfg, TFileDirectory &fs ):
        edm::BasicAnalyzer::BasicAnalyzer( cfg, fs ),
        src_( cfg.getParameter<edm::InputTag>( "src" ) ),
        genInfo_( cfg.getParameter<edm::InputTag>( "generatorInfo" ) ),
        processId_( cfg.getParameter<std::string>( "processId" ) ),
        lumiWeight_( cfg.getParameter<double>( "lumiWeight" ) ),
        maxCandPerEvent_( cfg.getParameter<int>( "maxCandPerEvent" ) ),
        sqrtS_( cfg.getUntrackedParameter<double>( "sqrtS", 13. ) ),
        intLumi_( cfg.getUntrackedParameter<double>( "intLumi",1000. ) ),
        nameTemplate_( cfg.getUntrackedParameter<std::string>( "nameTemplate", "$COLLECTION" ) ),
        dumpTrees_( cfg.getUntrackedParameter<bool>( "dumpTrees", false ) ),
        dumpWorkspace_( cfg.getUntrackedParameter<bool>( "dumpWorkspace", false ) ),
        workspaceName_( cfg.getUntrackedParameter<std::string>( "workspaceName", src_.label() ) ),
        dumpHistos_( cfg.getUntrackedParameter<bool>( "dumpHistos", false ) ),
        dumpGlobalVariables_( cfg.getUntrackedParameter<bool>( "dumpGlobalVariables", true ) ),
        classifier_( cfg.getParameter<edm::ParameterSet>( "classifierCfg" ) ),
        throwOnUnclassified_( cfg.exists("throwOnUnclassified") ? cfg.getParameter<bool>("throwOnUnclassified") : true ),
        globalVarsDumper_( 0 )        

    {
        if( cfg.getUntrackedParameter<bool>( "quietRooFit", false ) ) {
            RooMsgService::instance().setGlobalKillBelow( RooFit::WARNING );
        }
	    
        
        nPdfWeights_=0;
        dumpPdfWeights_=false;

        std::map<std::string, std::string> replacements;
        replacements.insert( std::make_pair( "$COLLECTION", src_.label() ) );
        replacements.insert( std::make_pair( "$PROCESS", processId_ ) );
        replacements.insert( std::make_pair( "$SQRTS", Form( "%1.0fTeV", sqrtS_ ) ) );
        nameTemplate_ = formatString( nameTemplate_, replacements );

        if( dumpGlobalVariables_ ) {
            globalVarsDumper_ = new GlobalVariablesDumper( cfg.getParameter<edm::ParameterSet>( "globalVariables" ) );
        }
       
        pdfWeightHistosBooked_=false;

        auto categories = cfg.getParameter<std::vector<edm::ParameterSet> >( "categories" );
        for( auto &cat : categories ) {
            auto label   = cat.getParameter<std::string>( "label" );
            auto subcats = cat.getParameter<int>( "subcats" );
            
            //------->
            //if any of the categories need pdfWeight dumping, we want to have the info in the WS
            //the logic of this code is that dumpPdfWeights_ is false or nPdfWeights_ is 0
            //we check to see what the category says. If any of the categories have more than 0 pdfWeight
            //and want dumpPdfWeights_ true, then the if fails and we don't check: the value stays true or >0
            //this relies on the fact that we want nPdfWeights the same for all cats.
            if (nPdfWeights_ == 0) {
            nPdfWeights_ = cat.exists("nPdfWeights") ?  cat.getParameter<int>( "nPdfWeights" ) : 0;
            }
            if (dumpPdfWeights_ == false ) {
            dumpPdfWeights_ = cat.exists("dumpPdfWeights")? cat.getParameter<bool>( "dumpPdfWeights" ) : false;
            }
            std::string classname = ( cat.exists("className") ? cat.getParameter<std::string>( "className" ) : "" );
            //<------
            
            auto name = nameTemplate_;
            if( ! label.empty() ) {
                name = replaceString( name, "$LABEL", label );
            } else {
                name = replaceString( replaceString( replaceString( name, "_$LABEL", "" ), "$LABEL_", "" ), "$LABEL", "" );
            }
            if( ! classname.empty() ) {
                name = replaceString( name, "$CLASSNAME", classname );
            } else {
                name = replaceString( replaceString( replaceString( name, "_$CLASSNAME", "" ), "$CLASSNAME_", "" ), "$CLASSNAME", "" );
            }
            
            KeyT key = classname;
            if( ! label.empty() ) {
                if( ! key.empty() ) { key += ":"; } // FIXME: define ad-hoc method with dedicated + operator
                key += label;
            }
            
            hasSubcat_[key] = ( subcats > 0 );
            auto &dumpers = dumpers_[key];
            if( subcats == 0 ) {
                name = replaceString( replaceString( replaceString( name, "_$SUBCAT", "" ), "$SUBCAT_", "" ), "$SUBCAT", "" );
                dumpers.push_back( dumper_type( name, cat, globalVarsDumper_ ) );
            } else {
                for( int isub = 0; isub < subcats; ++isub ) {
                    auto subcatname = replaceString( name, "$SUBCAT", Form( "%d", isub ) );
                    dumpers.push_back( dumper_type( subcatname, cat, globalVarsDumper_ ) );
                }
            }
        }

        if(dumpPdfWeights_){
            //if (cfg.exists("flashggPDFWeightObject")){
                pdfWeightToken_= cfg.getUntrackedParameter<edm::InputTag>("flashggPDFWeightObject",edm::InputTag("flashggPDFWeightObject"));
           // }
        }

        workspaceName_ = formatString( workspaceName_, replacements );
        if( dumpWorkspace_ ) {
            ws_ = fs.make<RooWorkspace>( workspaceName_.c_str(), workspaceName_.c_str() );
            dynamic_cast<RooRealVar *>( ws_->factory( "weight[1.]" ) )->setConstant( false );
            if (dumpPdfWeights_){
                for( int j=0; j<nPdfWeights_;j++ ) {
                    dynamic_cast<RooRealVar *>( ws_->factory( Form("pdfWeight_%d[1.]",j)) )->setConstant( false );
                }
            }
            RooRealVar* intLumi = new RooRealVar("IntLumi","IntLumi",intLumi_);
            //workspace expects intlumi in /pb
            intLumi->setConstant(); 
            // don't want this param to float in the fits at any point
            ws_->import(*intLumi);
        } else {
            ws_ = 0;
        }
        for( auto &dumpers : dumpers_ ) {
            for( auto &dumper : dumpers.second ) {
                if( dumpWorkspace_ ) {
                    dumper.bookRooDataset( *ws_, "weight", replacements);
                }
                if( dumpTrees_ ) {
                    TFileDirectory dir = fs.mkdir( "trees" );
                    dumper.bookTree( dir, "weight", replacements );
                }
                if( dumpHistos_ ) {
                    TFileDirectory dir = fs.mkdir( "histograms" );
                    dumper.bookHistos( dir, replacements );
                }
            }
        }
    }

    //// template<class C, class T, class U>
    //// CollectionDumper<C,T,U>::~CollectionDumper()
    //// {
    //// }

    template<class C, class T, class U>
        void CollectionDumper<C, T, U>::beginJob()
        {
        }

    template<class C, class T, class U>
        void CollectionDumper<C, T, U>::endJob()
        {
        }

    template<class C, class T, class U>
        double CollectionDumper<C, T, U>::eventWeight( const edm::EventBase &event )
        {
            double weight = 1.;
            if( ! event.isRealData() ) {
                edm::Handle<GenEventInfoProduct> genInfo;
                event.getByLabel( genInfo_, genInfo );

                weight = lumiWeight_;

                if( genInfo.isValid() ) {
                    const auto &weights = genInfo->weights();
                    // FIMXE store alternative/all weight-sets
                    if( ! weights.empty() ) {
                        weight *= weights[0];
                    }
                }
            }
            return weight;
        }

    template<class C, class T, class U>

        vector<double> CollectionDumper<C, T, U>::pdfWeights( const edm::EventBase &event )
        {   
            vector<double> pdfWeights;
            edm::Handle<vector<flashgg::PDFWeightObject> > WeightHandle;
            event.getByLabel( pdfWeightToken_, WeightHandle );


            for( unsigned int weight_index = 0; weight_index < (*WeightHandle).size(); weight_index++ ){
                std::vector<float> uncompressed = (*WeightHandle)[weight_index].uncompress();
                for( unsigned int j=0; j<(*WeightHandle)[weight_index].pdf_weight_container.size();j++ ) {
                    pdfWeights.push_back(uncompressed[j]);
                }
            }


            return pdfWeights;
        }



    template<class C, class T, class U>
        void CollectionDumper<C, T, U>::analyze( const edm::EventBase &event )
        {
            edm::Handle<collection_type> collectionH;
            event.getByLabel( src_, collectionH );
            const auto &collection = *collectionH;
            weight_ = eventWeight( event );

            if( globalVarsDumper_ ) { globalVarsDumper_->fill( event ); }
            if( dumpPdfWeights_){
                pdfWeights_ =pdfWeights( event );
            }

            int nfilled = maxCandPerEvent_;

            for( auto &cand : collection ) {
                auto cat = classifier_( cand );
                auto which = dumpers_.find( cat.first );

                if( which != dumpers_.end() ) {
                    int isub = ( hasSubcat_[cat.first] ? cat.second : 0 );
                    // FIXME per-candidate weights
                    which->second[isub].fill( cand, weight_, pdfWeights_, maxCandPerEvent_ - nfilled );
                    --nfilled;
                } else if( throwOnUnclassified_ ) {
                    throw cms::Exception( "Runtime error" ) << "could not find dumper for category [" << cat.first << "," << cat.second << "]"
                        << "If you want to allow this (eg because you don't want to dump some of the candidates in the collection)\n"
                        << "please set throwOnUnclassified in the dumper configuration\n";
                }
                if( ( maxCandPerEvent_ > 0 )  && nfilled == 0 ) { break; }
            }
        }

}

#endif // flashgg_CollectionDumper_h
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

