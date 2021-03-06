//
// Created by simonepanzeri on 24/11/2021.
//

#ifndef DEV_FDAPDE_DE_DATA_IMP_H
#define DEV_FDAPDE_DE_DATA_IMP_H

template<UInt ndim>
DEData<ndim>::DEData(const std::vector<Point<ndim>>& data, const UInt& order, const VectorXr& fvec, Real heatStep,
                     UInt heatIter, const std::vector<Real>& lambda, const UInt& nfolds, const UInt& nsim,
                     const std::vector<Real>& stepProposals, Real tol1, Real tol2, bool print, UInt search) :
        data_(data), order_(order), fvec_(fvec), heatStep_(heatStep), heatIter_(heatIter), lambda_(lambda),
        Nfolds_(nfolds), nsim_(nsim), stepProposals_(stepProposals), tol1_(tol1), tol2_(tol2), print_(print),
        search_(search) {}

template<UInt ndim>
void DEData<ndim>::printData(std::ostream & out) const
{
    for(int i = 0; i < data_.size(); i++)
    {
        out << data_[i] << std::endl;
    }
}

DEData_time::DEData_time(const std::vector<Real>& data_time, const std::vector<Real>& lambda_time) :
    data_time_(data_time), lambda_time_(lambda_time) {}

/*
template<UInt ndim>
void DEData_time::createMap(const std::vector<Point<ndim>>& data) {
    static constexpr Real eps = std::numeric_limits<Real>::epsilon(), tolerance = 100 * eps;

    data_time_noD_.clear();
    data_noD_.clear();

    //extraction of non duplicated times
    std::set<Real> set_time_noD(data_time_.cbegin(),data_time_.cend());
    data_time_noD_.resize(set_time_noD.size());
    std::copy(set_time_noD.begin(), set_time_noD.end(), data_time_noD_.begin());
    //std::vector<Real> v(set_time_noD.begin(), set_time_noD.end());
    //data_time_noD_.swap(v);
    set_time_noD.clear();

    //creation of the map which contains as keys the IDs of different spatial points and as values the indices of
    //the time instants point in data_time_noD_ (this structure is needed to build one of the FEmatrices)
    std::set<UInt> set_helper;
    //for(size_t i = data.front().id(); i <= data.back().id(); ++i) { // progressive ids (not considering out-of-range ids)
    UInt i;
    for(size_t n = 0; n < data.size(); ++n) {
        i = data[n].id();
        if(!isAlready(i, set_helper)) {insert(i, data_time_[n]); set_helper.insert(i);}
        for(size_t j = n + 1; j < data.size(); ++j){
            if (data[n].dist(data[j]) < tolerance and !isAlready(data[j].id(), set_helper))
            {insert(i, data_time_[j]); set_helper.insert(data[j].id());}
        }
    }
    set_helper.clear();
}

std::vector<UInt> DEData_time::getID_noD() const {
    std::vector<UInt> noDup;
    noDup.reserve(data_noD_.size());
    for(const auto &s : data_noD_)
        noDup.push_back(s.first);
    return noDup;
}

template <UInt ndim>
const std::set<UInt>& DEData_time<ndim>::getTimesIndices(UInt id) const {
    return data_noD_.find(id)->second;
};
*/

void DEData_time::setTimes2Locations() {
    static constexpr Real eps = std::numeric_limits<Real>::epsilon(), tolerance = 100 * eps;

    //extraction of non duplicated times
    std::set<Real> set_times(data_time_.cbegin(),data_time_.cend());
    times_.resize(set_times.size());
    (void)std::copy(set_times.begin(), set_times.end(), times_.begin());
    //std::vector<Real> v(set_time_noD.begin(), set_time_noD.end());
    //data_time_noD_.swap(v);
    set_times.clear();

    //creation of the data structure containing for each time index i, the location indices that are observed at time i
    //(ordered as data appear in the clean dataset).
    Times2Locations_.resize(times_.size());
    for(UInt t = 0; t < times_.size(); ++t)
    {
        for (UInt d = 0; d < data_time_.size(); ++d) {
            if (data_time_[d] == times_[t])
                Times2Locations_[t].push_back(d);
        }
    }
}

void DEData_time::printData(std::ostream& out) const {
    for(int i = 0; i < data_time_.size(); i++)
    {
        out << data_time_[i] << std::endl;
    }
}

void DEData_time::printTimes2Locations(std::ostream& out) const {
    for(std::vector<std::vector<UInt>>::const_iterator it = Times2Locations_.cbegin(); it != Times2Locations_.cend(); ++it) {
        out << "time index: " << it-Times2Locations_.cbegin() << '\t' << "[location index]:";
        for(UInt i : *it)
            out << " " << i;
        out << std::endl;
    }
}

#endif //DEV_FDAPDE_DE_DATA_IMP_H
