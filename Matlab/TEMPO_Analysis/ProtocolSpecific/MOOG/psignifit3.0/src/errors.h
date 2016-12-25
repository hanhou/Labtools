/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#ifndef ERRORS_H
#define ERRORS_H

/** \brief base class for all errors */
class PsiError {
    public:
        const char *message;
        PsiError() : message("Unspecified PsiError") {}
        PsiError(const char* message) : message(message) {}
        ~PsiError() {}
};

/** \brief Error class for errors that arise if some feature is not implemented
 *
 * This error should also be raised if a feature is not meant to be implemented as for purely
 * virtual methods
 */
class NotImplementedError : public PsiError {
};

/** \brief Error class for errors that are related to unsuitable arguments to a method or function */
class BadArgumentError : public PsiError {
    public:
        BadArgumentError() {}
        BadArgumentError(const char* message) : PsiError(message) {}
};

/** \brief Error class for errors that are related to small or large indices (out of range) */
class BadIndexError : public PsiError {
};

#endif
