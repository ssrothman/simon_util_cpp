#ifndef SROTHMAN_SIMONTOOLS_CA_H
#define SROTHMAN_SIMONTOOLS_CA_H

#include "particlePair.h"

namespace simon{
    /*
     * if T has members other than pT, eta, phi
     * these will have to be handled separately
     */
    template <typename T>
    inline void merge_particles(const T& p1, const T& p2,
                                T& result){
        result.pt = p1.pt + p2.pt;
        result.eta = (p1.pt*p1.eta + p2.pt*p2.eta) / result.pt;
        result.phi = (p1.pt*p1.phi + p2.pt*p2.phi) / result.pt;
    }

    template <bool distances_squared>
    struct indexedPair {
        const simon::pairentry_resolved<distances_squared> * const entry;
        const unsigned i; //index of particle 1
        const unsigned j; //index of particle 2
        const unsigned bkp; //bookkeeping index
        indexedPair(const simon::pairentry_resolved<distances_squared> * const entry,
                     const unsigned i, const unsigned j,
                     const unsigned bkp) noexcept :
            entry(entry), i(i), j(j), bkp(bkp) {}

        indexedPair() noexcept :
            entry(nullptr), i(0), j(0), bkp(0) {}
    };

    template <typename T, size_t N, bool distances_squared>
    inline void runCA(
            const std::array<const T*, N>& particles,
            const std::array<indexedPair<distances_squared>, N*(N-1)/2>& dRs,
            T& merged, unsigned &bkpindex){

        const indexedPair<distances_squared>& minDR = *std::min_element(
            dRs.begin(), dRs.end(),
            [](const indexedPair<distances_squared>& a, const indexedPair<distances_squared>& b){
                return a.entry->floatDR < b.entry->floatDR;
            });

        merge_particles(*particles[minDR.i], *particles[minDR.j], merged);
        bkpindex = minDR.bkp;
    }

    /*
     * There is only one topology
     *               ----p1
     *              /
     *         ----< r
     *        /     \
     *       /       ----p2
     *  ----< R
     *       \
     *        \
     *         ----------p3
     *
     * This can be fully described by:
     * particles: {p1, p2, p3, r, R}
     * pairs: {p1p2, rp3}
     *
     * The caller already has
     * particles: {p1, p2, p3}
     * pairs: {p1p2}
     *
     * so we only need to return:
     * particles: {r, R}
     * pairs: {rp3}
     * as well as an index to indicate which particle is p3
     *
     * So in summary, the function signature is:
     * CA3(particles, pairs, r, R, rp3, bkpindex);
     *
     * INPUTS: 
     *  particles: {p1, p2, p3}
     *  pairs: {p1p2, p1p3, p2p3}
     *  NB these particle indices need not 
     *          align with the diagram
     *
     * OUTPUTS:
     *  particles: {r, R}
     *  p1p2: resolved pair entry for p1p2
     *  bkpindex: index of p3 in the input array
     *      such that the topology diagram is correct
     */
    template <typename T, bool distances_squared>
    inline void CA3(
            const std::array<const T*, 3>& particles,
            const simon::pairentry_resolved<distances_squared>& p1p2,
            const simon::pairentry_resolved<distances_squared>& p1p3,
            const simon::pairentry_resolved<distances_squared>& p2p3,
            T& r, T& R, 
            simon::pairentry_resolved<distances_squared>& rp3,
            unsigned& bkpindex){

        //use the bookkeeping index to track which
        //particle is NOT in the pair
        const std::array<indexedPair<distances_squared>, 3> indexed_dRs = {
            indexedPair<distances_squared>(&p1p2, 0, 1, 2),
            indexedPair<distances_squared>(&p1p3, 0, 2, 1),
            indexedPair<distances_squared>(&p2p3, 1, 2, 0)
        };

        runCA(particles, indexed_dRs, r, bkpindex);

        const T* p3 = particles[bkpindex];
        rp3 = simon::pairentry_resolved<distances_squared>(r, *p3);
        merge_particles(r, *p3, R);
    }

    /*
     * four-point topology is either CHAIN or SYMMETRIC
     *
     * SYMMETRIC looks like:
     *
     *               ----p1
     *              /
     *         ----< r2
     *        /     \
     *       /       ----p2
     *  ----< R
     *       \       ----p3
     *        \     /
     *         ----< r1
     *              \
     *               ----p4
     *
     * The full information to describe this topology is
     * particles: {p1, p2, p3, p4, r1, r2, R}
     * pairs: {p1p2, p3p4, r1r2}
     *
     * NB the order of (p1, p2) and (p3, p4) is arbitrary
     *
     * CHAIN looks like:
     *
     *                       ----p1
     *                      /
     *                 ----< r2
     *                /     \
     *               /       ----p2
     *          ----< r1
     *         /     \
     *        /       -----------p3
     *   ----< R
     *        \
     *         \
     *          -----------------p4
     *
     * The full information to describe this topology is
     * particles: {p1, p2, p3, p4, r1, r2, R}
     * pairs: {p1p2, r2p3, r1p4}
     *
     * The caller already has
     * particles: {p1, p2, p3, p4}
     * pairs: {p1p2, p1p3, p1p4, p2p3, p2p4, p3p4}
     *
     * NB the order of (p1, p2) is arbitrary, 
     * but (p3, p4) is not
     *
     * so we only need to return:
     * for SYMMETRIC:
     *  particles: {r1, r2, R}
     *  pairs: {r1r2}
     *  four indices to indicate the order of p1, p2, p3, p4
     *      in the input array
     *
     * for CHAIN:
     *  particles: {r1, r2, R}
     *  pairs: {r1p4, r2p3}
     *  four indices to indicate the order of p1, p2, p3, p4
     *      in the input array
     *
     *
     * it's kinda awkward that there are two 
     * different return signatures, expecially
     * because part of the task of this function
     * is to determine which topology we have
     *
     * I think the best function signature is:
     * CA(particles, pairs, r1, r2, R, {return pairs}, topologyflag, {return indices});
     *
     * INPUTS:
     *  particles: {p1, p2, p3, p4}
     *  p1p2, p1p3, p1p4, p2p3, p2p4, p3p4: pair entries
     *  NB these particle indices need not
     *        align with the diagram
     *
     * OUTPUTS:
     *  r1, r2, R: merged particles
     *  return_pairs: {r1r2} or {r2p3, r1p4}
     *      depending on the topology
     *  topologyflag: CHAIN or SYMMETRIC
     *  indices: indices into input particles array
     *      such that p1, p2, p3, p4 are in the correct order
     *      as seen in the diagram above
     *  return_orig_pairs: pointers to {p1p2, p3p4} or {p1p2}
     *      depending on the topology
     *      these are just pointers to the input objects
     *      but it is useful to keep track of them
     *      and know which ones they are
     *      out of the six inputs 
     */
    constexpr int CHAIN = 0;
    constexpr int SYMMETRIC = 1;

    template <typename T, bool distances_squared>
    inline void CA4(
            const std::array<const T*, 4>& particles,
            const simon::pairentry_resolved<distances_squared>& p1p2,
            const simon::pairentry_resolved<distances_squared>& p1p3,
            const simon::pairentry_resolved<distances_squared>& p1p4,
            const simon::pairentry_resolved<distances_squared>& p2p3,
            const simon::pairentry_resolved<distances_squared>& p2p4,
            const simon::pairentry_resolved<distances_squared>& p3p4,
            T& r1, T& r2, T& R,
            std::vector<simon::pairentry_resolved<distances_squared>>& return_pairs,
            int& topologyflag,
            std::array<unsigned, 4>& return_indices,
            std::vector<const simon::pairentry_resolved<distances_squared>*>& return_orig_pairs){
        
        //use the bookkeeping index to track the index
        //of the pair that does not overlap 
        //with the current pair
        const std::array<indexedPair<distances_squared>, 6> indexed_dRs = {
            indexedPair<distances_squared>(&p1p2, 0, 1, 5),
            indexedPair<distances_squared>(&p1p3, 0, 2, 4),
            indexedPair<distances_squared>(&p1p4, 0, 3, 3),
            indexedPair<distances_squared>(&p2p3, 1, 2, 2),
            indexedPair<distances_squared>(&p2p4, 1, 3, 1),
            indexedPair<distances_squared>(&p3p4, 2, 3, 0)
        };

        unsigned bkp;
        runCA(particles, indexed_dRs, r2, bkp);

        //bkp finds us the opposite pair
        //and therefore which particles are p3 and p4
        const indexedPair<distances_squared>& opposite = indexed_dRs[bkp];
        const T* const p3 = particles[opposite.i];
        const T* const p4 = particles[opposite.j];
        return_indices[2] = opposite.i;
        return_indices[3] = opposite.j;
        //thanks to how nicely the bookkeeping index lines up
        //we actually also know the index of the merged pair
        const indexedPair<distances_squared>& merged = indexed_dRs[5-bkp];
        return_indices[0] = merged.i;
        return_indices[1] = merged.j;

        //the first merged pair is p1p2,
        //so we also want to store than in return_orig_pairs
        return_orig_pairs.resize(1);
        return_orig_pairs[0] = merged.entry;

        //now we set up to call CA3
        const std::array<const T*, 3> three_particles = {
            p3, p4, &r2
        };
        simon::pairentry_resolved<distances_squared> r2p3(r2, *p3);
        simon::pairentry_resolved<distances_squared> r2p4(r2, *p4);

        /*
         * CA3 will give us either r1r2 or r1p4
         * depending on whether we are chain or symmetric
         * in either case the pair wants to be the first one
         * in the return_pairs vector
         *
         * either way r1 and R will be correct
         * and we should be able to use the bookkeeping index
         * to determine which topology we have
         */
        return_pairs.resize(1);
        CA3(
            three_particles,
            *opposite.entry, r2p3, r2p4,
            r1, R, 
            return_pairs[0],
            bkp
        );

        if (bkp == 2){
            topologyflag = SYMMETRIC;
            //then we also care about which is p2p3 
            //for the return_orig_pairs
            return_orig_pairs.push_back(opposite.entry);
        } else {
            topologyflag = CHAIN;
            //then we need to return an extra pair
            return_pairs.push_back(r2p3);
            //and also the order of p3, p4 matters
            if (bkp == 0){
                std::swap(return_indices[2], return_indices[3]);
            }
        }
    }
};

#endif
