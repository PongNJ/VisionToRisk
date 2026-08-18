require("dotenv").config({ path: "google_key.env" });

function loadGoogleApiKeys() {
  return Object.entries(process.env)
    .filter(([name, value]) => /^GOOGLE_API_KEY_\d+$/.test(name) && value)
    .sort(([a], [b]) => a.localeCompare(b, undefined, { numeric: true }))
    .map(([, value]) => value);
}

function chooseApiKey() {
  const apiKeys = loadGoogleApiKeys();
  if (apiKeys.length === 0) {
    throw new Error("No GOOGLE_API_KEY_* environment variables are configured.");
  }
  return apiKeys[Math.floor(Math.random() * apiKeys.length)];
}

module.exports = {
  chooseApiKey,
  loadGoogleApiKeys,
};
