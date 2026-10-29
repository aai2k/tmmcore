/**
 * The names of 0.4, kept through 0.5.x with their old arguments and results,
 * and removed in 0.6.0. The first call of each prints one warning, naming the
 * call that replaces it when there is one.
 */

const warned = new Set();

/** Prints the warning for `name` the first time it is called, and never again. */
export function warnOnce(name, replacement = null) {
    if (warned.has(name)) return;
    warned.add(name);
    const instead = replacement ? `; use ${replacement}` : '';
    console.warn(`tmmcore: ${name} is deprecated and will be removed in 0.6.0${instead}`);
}

/**
 * `fn` under its old name: the same call, with the warning on first use, and
 * the name and parameter count the function had in 0.4.
 */
export function deprecated(name, replacement, fn) {
    const wrapper = function (...args) {
        warnOnce(name, replacement);
        return fn.apply(this, args);
    };
    return Object.defineProperties(wrapper, { name: { value: name }, length: { value: fn.length } });
}
